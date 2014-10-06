using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SequenceTag;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.Backend.Data.Spectrometry
{

    public class ChargeLcScanCluster
    {
        public double Score { get; private set; }

        public List<ChargeLcScanCell> Members;
        public double[] ObservedEnvelope;

        public double HighestIntensity
        {
            get 
            {
                return Members.Count < 1 ? 0.0 : Members[0].Intensity;
            }
        }
        

        public ChargeLcScanCluster()
        {
            Score = -1;
            Members = new List<ChargeLcScanCell>();
        }

        public void Add(ChargeLcScanCell cell, double[,,] data, double[] envelopePdf, double _score = -1)
        {
            var nEnvelope = data.GetLength(2);

            if (ObservedEnvelope == null)
                ObservedEnvelope = new double[nEnvelope];

            Members.Add(cell);

            for (var k = 0; k < nEnvelope; k++)
                ObservedEnvelope[k] += data[cell.Row, cell.Col, k];

            Score = _score >= 0 ? _score : ComputeScore(envelopePdf, ObservedEnvelope);
        }

        public double GetNewScore(ChargeLcScanCell cell, double[,,] data, double[] envelopePdf)
        {
            var nEnvelope = data.GetLength(2);
            var tempEnvelope = new double[nEnvelope];

            for (var k = 0; k < nEnvelope; k++)
                tempEnvelope[k] = ObservedEnvelope[k] + data[cell.Row, cell.Col, k];

            var newScore = ComputeScore(envelopePdf, tempEnvelope);

            return newScore;
        }

        public static double ComputeScore(double[] envelopPdf, double[] observedEnvelope)
        {
            const double minIntensity = 1E-6;
            var normalizedEnv = new double[observedEnvelope.Length];

            var sum = observedEnvelope.Sum();
            bool foundZero = false;
            for (var k = 0; k < envelopPdf.Length; k++)
            {
                normalizedEnv[k] = observedEnvelope[k] / sum;
                if (normalizedEnv[k] < minIntensity)
                {
                    normalizedEnv[k] = minIntensity;
                    foundZero = true;
                }
            }

            if (!foundZero) return SimpleMath.GetKLDivergence(envelopPdf, normalizedEnv);

            sum = normalizedEnv.Sum();
            for (var k = 0; k < envelopPdf.Length; k++)
                normalizedEnv[k] = normalizedEnv[k] / sum;

            return SimpleMath.GetKLDivergence(envelopPdf, normalizedEnv);
        }
        
    }

    public class ChargeLcScanCell //: IComparable<ChargeLcScanCell>
    {
        public int Row;
        public int Col;
        public double Intensity;

        public ChargeLcScanCell(int row, int col, double intensity)
        {
            Row = row;
            Col = col;
            Intensity = intensity;
        }

        public int ChargeIndex
        {
            get { return Row; }
        }

        public int ScanIndex
        {
            get { return Col;  }
        }

        //public int CompareTo(ChargeLcScanCell other)
        //{
        //    return Intensity.CompareTo(other.Intensity);
        //}
    }

    public class ChargeLcScanMatrix
    {
        private LcMsRun _run;
        public int NumberOfScans { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }

        public double[,,] Data { get; private set; }
        public double[,] SumOverIsotope { get; private set; }

        public double[] EnvelopePdf;

        public int NRows
        {
            get { return MaxCharge - MinCharge + 1; }
        }

        public int NColumns
        {
            get { return NumberOfScans; }
        }
        
        public ChargeLcScanMatrix(LcMsRun run)
        {
            _run = run;
            MinCharge = 2;
            MaxCharge = 60;
            NumberOfScans = _run.GetScanNumbers(1).Count;
        }

        public void SetChargeRange(int minCharge, int maxCharge)
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
        }

        private void BuildMatrix(double proteinMass, int mzBinBitSize = 27)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            EnvelopePdf = new double[isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                EnvelopePdf[isotopeIndex] = isoEnv.Envolope[isotopeIndex]/isoEnv.Envolope.Sum();

            var comparer = new MzComparerWithBinning(mzBinBitSize);
            SumOverIsotope = new double[NRows, NColumns];
            Data = new double[NRows, NColumns, isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
            {
                for (var charge = MinCharge; charge <= MaxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var binNum = comparer.GetBinNumber(abundantIsotopeMz);

                    var mz = comparer.GetMz(binNum);
                    var mzNext = comparer.GetMz(binNum + 1);

                    var xic = _run.GetFullPrecursorIonExtractedIonChromatogram(mz, mzNext);

                    var intensities = xic.Select(p => p.Intensity);
                    var scanIndex = 0;
                    foreach (var intensity in intensities)
                    {
                        SumOverIsotope[charge - MinCharge, scanIndex] += intensity;
                        Data[charge - MinCharge, scanIndex, isotopeIndex] = intensity;
                        scanIndex++;
                    }
                }
            }
        }


        private double GetIntensityThreshold(double[,] matrix)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            for (var i = 0; i < NRows; i++)
                for (var j = 0; j < NColumns; j++)
                    if (matrix[i, j] > 0) nonZeroIntensities.Add(matrix[i, j]);

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[nonZeroIntensities.Count / 2];

            return intensityThreshold;
        }


        public IEnumerable<ChargeLcScanCluster> GetProbableChargeScanRegions(double proteinMass)
        {
            BuildMatrix(proteinMass);
            var clusters = Segmentation(SumOverIsotope);

            return clusters;
        }
        
        public List<ChargeLcScanCluster> Segmentation(double[,] matrix)
        {
            var intensityThreshold = GetIntensityThreshold(matrix);
            var tempMap = new byte[NRows, NColumns];

            var seedCells = new List<ChargeLcScanCell>();

            
            for(var i = 0; i < NRows; i++)
                for (var j = 0; j < NColumns; j++)
                    if (matrix[i, j] > intensityThreshold)
                    {
                        tempMap[i, j] = 1;
                        seedCells.Add(new ChargeLcScanCell(i, j, matrix[i, j]));
                    }

            //seedCells.Sort();

            var orderedSeedCells = seedCells.OrderByDescending(cell => cell.Intensity);

            //1) Label start as 1
	        //2) Get neighbor cells based on start point ->insert into a list L
	        //3) Label neighbor as 2
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seed in orderedSeedCells)
            {
                if (tempMap[seed.Row, seed.Col] != 1) continue;

                var newCluster = new ChargeLcScanCluster();
                var neighbors = new Queue<ChargeLcScanCell>();

                // pick a seed cell (i, j)
                neighbors.Enqueue(seed);
                
                newCluster.Add(seed, Data, EnvelopePdf);
                tempMap[seed.Row, seed.Col] = 2;  

                do
                {
                    var cell = neighbors.Dequeue();
                    /*
                    newCluster.Add(cell, Data, EnvelopePdf);
                    tempMap[cell.Row, cell.Col] = 2;  
                    */
                    //collect neighbor cells around the focused pixel
                    for (var k = cell.Row - 1; k < cell.Row + 2; k++)
                    {
                        for (var l = cell.Col - 1; l < cell.Col + 2; l++)
                        {
                            if (k < 0 || l < 0 || k >= NRows || l >= NColumns || tempMap[k, l] == 2 || matrix[k, l] < newCluster.HighestIntensity * 0.1) continue;

                            var newMember = new ChargeLcScanCell(k, l, matrix[k, l]);
                            var newScore = newCluster.GetNewScore(newMember, Data, EnvelopePdf);

                            if (newScore >= newCluster.Score)
                            {
                                neighbors.Enqueue(newMember);
                                newCluster.Add(newMember, Data, EnvelopePdf, newScore);
                                tempMap[k, l] = 2;
                            }
                        }
                    }
                }
                while (neighbors.Count > 0);

                //if (newCluster.Score == Double.PositiveInfinity) continue;

                clusters.Add(newCluster);                 
            }

            return clusters.OrderBy(x => x.Score).ToList();
        }

        /*
        private ChargeLcScanCluster CreateCluster(List<ChargeLcScanCell> cells)
        {
            var startScanNum = NumberOfScans;
            var endScanNum = 0;

            var maxC = 0;
            var minC = MaxCharge - MinCharge;

            foreach (var cell in cells)
            {
                if (cell.ScanIndex > endScanNum) endScanNum = cell.ScanIndex;
                if (cell.ScanIndex < startScanNum) startScanNum = cell.ScanIndex;

                if (cell.ChargeIndex > maxC) maxC = cell.ChargeIndex;
                if (cell.ChargeIndex < minC) minC = cell.ChargeIndex;
            }

            var cluster = new ChargeLcScanCluster(MinCharge, minC, maxC, startScanNum, endScanNum);
            
            cluster.Score1 = ComputeScore1(cluster);
            cluster.Score2 = ComputeScore2(cluster);
            cluster.Score3 = ComputeScore3(cluster);

            cluster.Score = cluster.Score1 + cluster.Score2 + cluster.Score3;
            
            return cluster;
        }
        
        private double ComputeScore1(ChargeLcScanCluster cluster)
        {
            if (NumberOfIsotopeEnvelop < 2) return 0.0;
            
            var ipOverIsotopes = new double[NumberOfIsotopeEnvelop];

            for (var k = 0; k < NumberOfIsotopeEnvelop; k++)
            {
                for (var i = 0; i < NRows; i++)
                    for (var j = cluster.MinScanIndex; j <= cluster.MaxScanIndex; j++)
                        ipOverIsotopes[k] += Data[i, j, k];
            }

            return SimpleMath.GetCorrelation(_isotopeEnvelop, ipOverIsotopes);
        }

        private double ComputeScore2(ChargeLcScanCluster cluster)
        {
            if (cluster.ChargeCount < 2 || cluster.ScanCount < 2) return 0.0;
            
            var ipOverChargeTime = new double[cluster.ChargeCount][];

            for (int i = 0; i < cluster.MaxChargeIndex - cluster.MinChargeIndex + 1; i++)
                ipOverChargeTime[i] = new double[cluster.MaxScanIndex - cluster.MinScanIndex + 1];
            
            for (var i = cluster.MinChargeIndex; i <= cluster.MaxChargeIndex; i++)
                for (var j = cluster.MinScanIndex; j <= cluster.MaxScanIndex; j++)
                    for (var k = 0; k < NumberOfIsotopeEnvelop; k++)
                        ipOverChargeTime[i - cluster.MinChargeIndex][j - cluster.MinScanIndex] += Data[i, j, k];

            // do we want to calculate correlation between the most two intense charges?
            //var correlations = new List<double>();
            //for (var i = 0; i < ipOverChargeTime.Length; i++)
            //    for (var j = i + 1; j < ipOverChargeTime.Length; j++)
            //        correlations.Add( SimpleMath.GetCorrelation(ipOverChargeTime[i], ipOverChargeTime[j]) );
            //return SimpleMath.GetSampleMean(correlations.ToArray());
            var temp = new List<double>();
            for (var i = 0; i < cluster.ChargeCount; i++)
                temp.Add(ipOverChargeTime[i].Sum());

            var tempSort = temp.Select((x, i) => new KeyValuePair<double, int>(x, i)).OrderByDescending(x => x.Key).ToList();

            return SimpleMath.GetCorrelation(ipOverChargeTime[tempSort[0].Value], ipOverChargeTime[tempSort[1].Value]);
        }

        private double ComputeScore3(ChargeLcScanCluster cluster)
        {
            if (NumberOfIsotopeEnvelop < 2 || cluster.ScanCount < 2) return 0.0;
            
            var ipOverisotopeTime = new double[NumberOfIsotopeEnvelop][];

            for (int i = 0; i < NumberOfIsotopeEnvelop; i++)
                ipOverisotopeTime[i] = new double[cluster.ScanCount];

            for (var k = 0; k < NumberOfIsotopeEnvelop; k++)
                for (var j = cluster.MinScanIndex; j <= cluster.MaxScanIndex; j++)
                    for (var i = cluster.MinChargeIndex; i <= cluster.MaxChargeIndex; i++)
                        ipOverisotopeTime[k][j - cluster.MinScanIndex] += Data[i, j, k];
            
            //var correlations = new List<double>();
            //for (var i = 0; i < ipOverisotopeTime.Length; i++)
            //    for (var j = i+1; j < ipOverisotopeTime.Length; j++)
            //        correlations.Add(SimpleMath.GetCorrelation(ipOverisotopeTime[i], ipOverisotopeTime[j]));
            //return SimpleMath.GetSampleMean(correlations.ToArray());

            return SimpleMath.GetCorrelation(ipOverisotopeTime[0], ipOverisotopeTime[1]);
        }
        */
        /*
        public IEnumerable<ChargeLcScanCluster> GetProbableChargeScanRegions()
        {
            var clusters = Segmentation();

            for (var i = 0; i < clusters.Count; i++)
            {
                var largestSegment = clusters[i].Value;

                var startScanNum = NumberOfLcScans;
                var endScanNum = 0;

                var maxC = MinCharge;
                var minC = MaxCharge;

                foreach (var cell in largestSegment)
                {
                    if (cell.ScanIndex > endScanNum) endScanNum = cell.ScanIndex;
                    if (cell.ScanIndex < startScanNum) startScanNum = cell.ScanIndex;

                    if (cell.ChargeIndex + MinCharge > maxC) maxC = cell.ChargeIndex + MinCharge;
                    if (cell.ChargeIndex + MinCharge < minC) minC = cell.ChargeIndex + MinCharge;
                }

                if (i > 10 || (i > 1 && clusters[i].Key < 10)) yield break;

                yield  return new ChargeLcScanCluster(minC, maxC, startScanNum, endScanNum);
            }
        }*/
    }
}
