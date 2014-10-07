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
    public class ChargeLcScanMatrix : ILcMsMap
    {
        private LcMsRun _run;
        private readonly int _numberOfScans;
        private int _minCharge;
        private int _maxCharge;
        private double[,,] _data;
        private double[,] _sumOverIsotope;
        private double[] _envelopePdf;

        private const double ScoreThreshold = 0.2;
        private readonly int[] _ms1ScanNums;
        int NRows { get { return _maxCharge - _minCharge + 1; } }
        int NColumns { get { return _numberOfScans; } }
        
        private readonly MzComparerWithBinning _comparer;
        private Dictionary<int, double[]> _cachedXicIntensity;

        public ChargeLcScanMatrix(LcMsRun run, int numBits = 27)
        {
            _minCharge = 2;
            _maxCharge = 60;
            _ms1ScanNums = run.GetMs1ScanVector();
            _numberOfScans = _ms1ScanNums.Length;
            _comparer = new MzComparerWithBinning(numBits);

            _cachedXicIntensity = new Dictionary<int, double[]>();
            _run = run;

            /*
            _startBinNum = _comparer.GetBinNumber(minMz);
            var endBinNum = _comparer.GetBinNumber(maxMz);
            _chachedXicIntensity = new double[endBinNum - _startBinNum + 1][];
            
            //pre-generation 
            for (var binNum = _startBinNum; binNum <= endBinNum; binNum++)
            {
                var mzStart = _comparer.GetMzStart(binNum);
                var mzEnd = _comparer.GetMzEnd(binNum);

                _chachedXicIntensity[binNum - _startBinNum] = run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);
            }*/
        }

        public void SetChargeRange(int minCharge, int maxCharge)
        {
            _minCharge = minCharge;
            _maxCharge = maxCharge;
        }

        public IEnumerable<List<Tuple<int, int>>> GetProbableChargeScanRegions(double proteinMass)
        {
            BuildMatrix(proteinMass);
            var clusters = Segmentation(_sumOverIsotope);
            var filtered = clusters.Where(cluster => cluster.GetScore() < ScoreThreshold).ToList();

            return filtered.Select(cluster => cluster.Members.Select(member => new Tuple<int, int>(member.Row + _minCharge, _ms1ScanNums[member.Col])).ToList()).ToList();
        }

        public IEnumerable<List<Tuple<int, int, double>>> GetAllChargeScanRegions(double proteinMass)
        {
            BuildMatrix(proteinMass);
            var clusters = Segmentation(_sumOverIsotope);
            return clusters.Select(cluster => cluster.Members.Select(member => new Tuple<int, int, double>(member.Row, member.Col, cluster.GetScore())).ToList()).ToList();
        }


        private void BuildMatrix(double proteinMass)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            _envelopePdf = new double[isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                _envelopePdf[isotopeIndex] = isoEnv.Envolope[isotopeIndex]/isoEnv.Envolope.Sum();
            
            _sumOverIsotope = new double[NRows, NColumns];
            _data = new double[NRows, NColumns, isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
            {
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var binNum = _comparer.GetBinNumber(abundantIsotopeMz);
                    double[] intensities = null;

                    if (_cachedXicIntensity.ContainsKey(binNum))
                    {
                        intensities = _cachedXicIntensity[binNum];
                    }
                    else
                    {
                        var mz = _comparer.GetMzAverage(binNum);
                        var mzNext = _comparer.GetMzAverage(binNum + 1);
                        var xic = _run.GetFullPrecursorIonExtractedIonChromatogram(mz, mzNext);
                        intensities = xic.Select(p => p.Intensity).ToArray();
                        _cachedXicIntensity.Add(binNum, intensities);
                    }

                    for(var scanIndex = 0; scanIndex < intensities.Length; scanIndex++)
                    {
                        _sumOverIsotope[charge - _minCharge, scanIndex] += intensities[scanIndex];
                        _data[charge - _minCharge, scanIndex, isotopeIndex] = intensities[scanIndex];
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
       
        
        private IEnumerable<ChargeLcScanCluster> Segmentation(double[,] matrix)
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
                
                newCluster.Add(seed, _data, _envelopePdf);
                tempMap[seed.Row, seed.Col] = 2;  

                do
                {
                    var cell = neighbors.Dequeue();
                    //collect neighbor cells around the focused pixel
                    for (var k = cell.Row - 1; k < cell.Row + 2; k++)
                    {
                        for (var l = cell.Col - 1; l < cell.Col + 2; l++)
                        {
                            if (k < 0 || l < 0 || k >= NRows || l >= NColumns || tempMap[k, l] == 2 || matrix[k, l] < newCluster.HighestIntensity * 0.1) continue;

                            var newMember = new ChargeLcScanCell(k, l, matrix[k, l]);
                            var newScore = newCluster.GetNewScore(newMember, _data, _envelopePdf);

                            if (newScore >= newCluster.GetScore())
                            {
                                neighbors.Enqueue(newMember);
                                newCluster.Add(newMember, _data, _envelopePdf, newScore);
                                tempMap[k, l] = 2;
                            }
                        }
                    }
                }
                while (neighbors.Count > 0);

                if (newCluster.Members.Count == 1)
                {
                    tempMap[newCluster.Members[0].Row, newCluster.Members[0].Col] = 1;
                }
                else clusters.Add(newCluster);                 
            }

            return clusters.OrderBy(x => x.GetScore()).ToList();
        }
        }

    class ChargeLcScanCluster
    {
        private double _score;

        internal readonly List<ChargeLcScanCell> Members;
        private double[] _observedEnvelope;

        internal ChargeLcScanCluster()
        {
            _score = -1;
            Members = new List<ChargeLcScanCell>();
        }

        internal double GetScore()
        {
            return _score;
        }

        internal void Add(ChargeLcScanCell cell, double[,,] data, double[] envelopePdf, double _score = -1)
        {
            var nEnvelope = data.GetLength(2);

            if (_observedEnvelope == null)
                _observedEnvelope = new double[nEnvelope];

            Members.Add(cell);

            for (var k = 0; k < nEnvelope; k++)
                _observedEnvelope[k] += data[cell.Row, cell.Col, k];

            this._score = _score >= 0 ? _score : ComputeScore(envelopePdf, _observedEnvelope);
        }

        internal double HighestIntensity
        {
            get { return Members.Count < 1 ? 0.0 : Members[0].Intensity; }
        }

        internal double GetNewScore(ChargeLcScanCell cell, double[,,] data, double[] envelopePdf)
        {
            var nEnvelope = data.GetLength(2);
            var tempEnvelope = new double[nEnvelope];

            for (var k = 0; k < nEnvelope; k++)
                tempEnvelope[k] = _observedEnvelope[k] + data[cell.Row, cell.Col, k];

            var newScore = ComputeScore(envelopePdf, tempEnvelope);

            return newScore;
        }

        internal double ComputeLinearCorrelation(double[] envelopPdf)
        {
            return SimpleMath.GetCorrelation(envelopPdf, _observedEnvelope);
        }

        internal double ComputeScore(double[] envelopPdf, double[] observedEnvelope)
        {
            //return 1 - SimpleMath.GetCorrelation(envelopPdf, _observedEnvelope);

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

    class ChargeLcScanCell
    {
        internal int Row;
        internal int Col;
        internal double Intensity;

        internal ChargeLcScanCell(int row, int col, double intensity)
        {
            Row = row;
            Col = col;
            Intensity = intensity;
        }
    }
}
