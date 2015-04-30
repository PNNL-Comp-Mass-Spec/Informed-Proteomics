using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using Constants = InformedProteomics.Backend.Data.Biology.Constants;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms1FeatureMatrix
    {
        public Ms1FeatureMatrix(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0, int numBits = 27)
        {
            Run = run;
            MinScanCharge  = minScanCharge;
            MaxScanCharge  = maxScanCharge;
            _maxThreadCount = maxThreadCount;
            Comparer        = new MzComparerWithBinning(numBits);
            _ms1PeakList    = new List<Ms1Peak>();
            Spectrums      = new List<Ms1Spectrum>();
            var ms1ScanNums = run.GetMs1ScanVector();
            NScans = ms1ScanNums.Length;

            for (var i = 0; i < Math.Min(ms1ScanNums.Length, ushort.MaxValue); i++)
            {
                var ms1Spec = run.GetMs1Spectrum(ms1ScanNums[i]);
                Spectrums.Add(ms1Spec);
                _ms1PeakList.AddRange((Ms1Peak[])ms1Spec.Peaks);
            }
            _ms1PeakList.Sort();

            CorrelationMap     = new double[MaxChargeLength][];
            DistanceMap        = new double[MaxChargeLength][];
            AccurateMass       = new double[MaxChargeLength][];
            FeatureMatrix      = new Ms1Peak[MaxChargeLength][][];
            _checkedOut         = new bool[MaxChargeLength][];

            for (var i = 0; i < MaxChargeLength; i++)
            {
                _checkedOut[i]          = new bool[NScans];
                CorrelationMap[i]      = new double[NScans];
                FeatureMatrix[i]       = new Ms1Peak[NScans][];
                DistanceMap[i]         = new double[NScans];
                AccurateMass[i]        = new double[NScans];

                for (var j = 0; j < NScans; j++)
                {
                    FeatureMatrix[i][j] = new Ms1Peak[MaxEnvelopeLength];
                }
            }
        }

        public IList<Ms1FeatureCluster> GetProbableClusters(int queryMassBinNum)
        {
            return GetProbableClusters(Comparer.GetMzAverage(queryMassBinNum));
        }

        public IList<Ms1FeatureCluster> GetProbableClusters(double queryMass)
        {
            SetQueryMass(queryMass);
            var clusters = FindClusters(true);
            return clusters;
        }

        public IList<Ms1FeatureCluster> GetAllClusters(double queryMass)
        {
            SetQueryMass(queryMass);
            var clusters = FindClusters(false);
            return clusters;
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(Ms1FeatureCluster cluster)
        {
            var ms2ScanNumSet = new List<int>();
            for (var charge = MinScanCharge; charge <= MaxScanCharge; charge++)
            {
                var mostAbuMz = Ion.GetIsotopeMz(cluster.RepresentativeMass, charge, TheoreticalEnvelope.GetMostAbundantIsotope().Index);
                var ms2ScanNums = Run.GetFragmentationSpectraScanNums(mostAbuMz);
                ms2ScanNumSet.AddRange(ms2ScanNums.Where(sn => sn > cluster.MinScanNum && sn < cluster.MaxScanNum));
            }

            return ms2ScanNumSet.Distinct();
        }

        public Ms1Spectrum GetSpectrum(int ms1ScanNum)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var idx = Array.BinarySearch(ms1ScanNums, ms1ScanNum);
            return idx < 0 ? null : Spectrums[idx];
        }

        public readonly MzComparerWithBinning Comparer;
        public readonly static Tolerance MzTolerance = new Tolerance(5);

        public LcMsRun Run;
        public List<Ms1Spectrum> Spectrums;
        public int MinScanCharge;
        public int MaxScanCharge;
        public int NScans;
        private readonly int _maxThreadCount;
        
        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;
        protected const int MaxChargeLength = 40;
        private readonly bool[][] _checkedOut;

        protected readonly double[][] CorrelationMap;
        protected readonly double[][] DistanceMap;
        protected readonly double[][] AccurateMass;
        protected readonly Ms1Peak[][][] FeatureMatrix;

        protected IsotopeList TheoreticalEnvelope;
        protected int[] Rows;
        protected int[] Cols;
        public IntRange ChargeRange;

        protected int MinSearchMassBin;
        protected int MaxSearchMassBin;
        protected double QueryMass;

        protected void SetQueryMass(double queryMass)
        {
            QueryMass = queryMass;
            ChargeRange = GetScanChargeRange(QueryMass);
            //var queryMassBinNum     = Comparer.GetBinNumber(QueryMass);
            TheoreticalEnvelope = new IsotopeList(queryMass, MaxEnvelopeLength);
        }

        protected IntRange GetScanChargeRange(double mass)
        {
            if (mass < 5000.0d) return new IntRange(MinScanCharge, 10);

            var chargeLb = (int)Math.Max(MinScanCharge, Math.Floor((13.0 / 2.5) * (mass / 10000d) - 0.6));
            var chargeUb = (int)Math.Min(MaxScanCharge, Math.Ceiling(18 * (mass / 10000d) + 8));

            if (chargeUb - chargeLb + 1 > MaxChargeLength) chargeUb = chargeLb + MaxChargeLength - 1;

            return new IntRange(chargeLb, chargeUb);
        }

        protected void BuildFeatureMatrix()
        {
            var queryMassBinNum     = Comparer.GetBinNumber(QueryMass);

            var options             = new ParallelOptions();
            var observedRows        = new bool[ChargeRange.Length];
            var observedCols        = new bool[NScans];
            var rows                = Enumerable.Range(0, ChargeRange.Length);

            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            var mostAbuIsotopeInternalIndex = TheoreticalEnvelope.SortedIndexByIntensity[0];
            var mostAbuIsotopeIndex = TheoreticalEnvelope[mostAbuIsotopeInternalIndex].Index;

            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;

            Parallel.ForEach(rows, options, row =>
            {
                Array.Clear(CorrelationMap[row], 0, NScans);
                Array.Clear(_checkedOut[row], 0, NScans);
                Array.Clear(DistanceMap[row], 0, NScans);
                Array.Clear(AccurateMass[row], 0, NScans);

                for (var col = 0; col < NScans; col++)
                {
                    Array.Clear(FeatureMatrix[row][col], 0, FeatureMatrix[row][col].Length);
                }

                var charge = row + ChargeRange.Min;
                for (var k = 0; k < TheoreticalEnvelope.Count; k++)
                {
                    var i = TheoreticalEnvelope.SortedIndexByIntensity[k]; // internal isotope index
                    var isotopeIndex = TheoreticalEnvelope[i].Index;

                    var isotopeMzLb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum - 1), charge, isotopeIndex);
                    var isotopeMzUb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum + 1), charge, isotopeIndex);
                    
                    //if (isotopeMzLb < _run.MinMs1Mz || isotopeMzUb > _run.MaxMs1Mz) continue;
                    if (isotopeMzLb < minMs1Mz || isotopeMzUb > maxMs1Mz) continue;
                    var st = _ms1PeakList.BinarySearch(new Ms1Peak(isotopeMzLb, 0, 0));

                    if (st < 0) st = ~st;

                    for (var j = st; j < _ms1PeakList.Count; j++)
                    {
                        var ms1Peak         = _ms1PeakList[j];
                        if (ms1Peak.Mz > isotopeMzUb) break;
                        var col = ms1Peak.Ms1SpecIndex;

                        if (k < 1) // if (k < 4)
                        {
                            if (FeatureMatrix[row][col][i] == null || ms1Peak.Intensity > FeatureMatrix[row][col][i].Intensity)
                            {
                                FeatureMatrix[row][col][i] = ms1Peak;
                                AccurateMass[row][col] = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                            }
                            // In case of top 3 isotopes, intense peaks are preferred
                            //if (_featureMatrix[row][col][i] == null || ms1Peak.Intensity > _featureMatrix[row][col][i].Intensity) _featureMatrix[row][col][i] = ms1Peak;

                            // accurate mass is derived from the most intense peak among top 3 isoptes
                            //if ((!(_accurateMass[row][col] > 0)) || _featureMatrix[row][col][i].Intensity > _featureMatrix[row][col][_highestPeakIsotopeInternalIndex[row][col]].Intensity)
                            //{
                                //_highestPeakIsotopeInternalIndex[row][col] = i;
                                //_accurateMass[row][col]         = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                            //}
                        }
                        else
                        {
                            var mostAbuPeak = FeatureMatrix[row][col][mostAbuIsotopeInternalIndex];
                            if (mostAbuPeak == null) continue;
                            var expectedPeakMz = mostAbuPeak.Mz + (Constants.C13MinusC12*(isotopeIndex - mostAbuIsotopeIndex))/charge;

                            //var mz = Ion.GetIsotopeMz(_accurateMass[row][col], charge, isotopeIndex);
                            if (Math.Abs(expectedPeakMz - ms1Peak.Mz) > MzTolerance.GetToleranceAsTh(ms1Peak.Mz)) continue;

                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (FeatureMatrix[row][col][i] != null)
                            {
                                var tmpPeak = FeatureMatrix[row][col][i];
                                var bc1 = FeatureMatrix[row][col].GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);
                                FeatureMatrix[row][col][i] = ms1Peak;
                                var bc2 = FeatureMatrix[row][col].GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);
                                if (bc1 < bc2) FeatureMatrix[row][col][i] = tmpPeak;
                            }
                            else
                            {
                                FeatureMatrix[row][col][i] = ms1Peak;
                            }
                        }
                    }
                }

                for (var col = 0; col < NScans; col++)
                {
                    if (!(AccurateMass[row][col] > 0)) continue;
                    var mostAbuPeakIntensity = FeatureMatrix[row][col][mostAbuIsotopeInternalIndex].Intensity;
                    var signalToNoiseRatio = mostAbuPeakIntensity / Spectrums[col].MedianIntensity;
                    if (signalToNoiseRatio > 1.4826 && FeatureMatrix[row][col].Count(p => p != null && p.Active) >= 3)
                    {
                        CorrelationMap[row][col] = FeatureMatrix[row][col].GetPearsonCorrelation(TheoreticalEnvelope.Envelope);
                        DistanceMap[row][col] = FeatureMatrix[row][col].GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);
                        if (!observedRows[row]) observedRows[row] = true;
                        if (!observedCols[col]) observedCols[col] = true;
                    }
                    else
                    {
                        AccurateMass[row][col] = 0d;
                        //_distanceMap[row][col] = 1.0d;
                    }
                }
            }// end or row for-loop
            );

            var temp = new List<int>();
            for (var i = 0; i < observedRows.Length; i++) if (observedRows[i]) temp.Add(i);
            Rows = temp.ToArray();
            temp.Clear();
            for (var i = 0; i < observedCols.Length; i++) if (observedCols[i]) temp.Add(i);
            Cols = temp.ToArray();

        }

        protected virtual IEnumerable<ObservedEnvelope> GetSeedCells()
        {
            var seedList = new List<KeyValuePair<double, ObservedEnvelope>>();
            foreach (var i in Rows)
            {
                foreach (var j in Cols)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;

                    var mostAbuPeakIntensity = FeatureMatrix[i][j][TheoreticalEnvelope.SortedIndexByIntensity[0]].Intensity;
                    var signalToNoiseRatio = mostAbuPeakIntensity / Spectrums[j].MedianIntensity;
                    if (signalToNoiseRatio < 3) continue;

                    if (FeatureMatrix[i][j].Count(p => p != null && p.Active) < 5) continue;

                    var bcDist = DistanceMap[i][j];
                    if (bcDist > 0.3 && CorrelationMap[i][j] < 0.5) continue;
                    var seedCell = new ObservedEnvelope(i, j, FeatureMatrix[i][j], TheoreticalEnvelope);
                    seedList.Add(new KeyValuePair<double, ObservedEnvelope>(bcDist, seedCell));
                }
            }
            return seedList.OrderBy(x => x.Key).Select(x => x.Value);
        }

        protected virtual double GetBcDistTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                if (normalizedElutionLen < 0.005) return 0.6;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.2;
                return 0.1;                
            }
            else if (QueryMass < 25000)
            {
                if (normalizedElutionLen < 0.005) return 1.0;
                if (normalizedElutionLen < 0.01) return 0.5;
                if (normalizedElutionLen < 0.02) return 0.25;
                return 0.2;                
            }
            else // > 25K 
            {
                if (normalizedElutionLen < 0.005) return 1.2;
                if (normalizedElutionLen < 0.01) return 0.8;
                if (normalizedElutionLen < 0.02) return 0.3;
                return 0.2;
            }
        }

        protected virtual double GetCorrTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                if (normalizedElutionLen < 0.005) return 0.3;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.6;
                return 0.7;
            }
            else if (QueryMass < 25000)
            {
                if (normalizedElutionLen < 0.005) return 0;
                if (normalizedElutionLen < 0.01) return 0.2;
                if (normalizedElutionLen < 0.02) return 0.4;
                return 0.6;
            }
            else // 25K
            {
                if (normalizedElutionLen < 0.005) return -1;
                if (normalizedElutionLen < 0.01) return 0.1;
                if (normalizedElutionLen < 0.02) return 0.4;
                return 0.5;
            }
        }

        protected virtual int GetScatteredChargeLength(int charge)
        {
            var scatteredLength = 1;
            if (charge > 40) scatteredLength = 12;
            else if (charge > 30) scatteredLength = 8;
            else if (charge > 20) scatteredLength = 4;
            else if (charge > 10) scatteredLength = 2;
            return scatteredLength;
        }

        private List<Ms1FeatureCluster> FindClusters(bool filtering = false)
        {
            BuildFeatureMatrix(); // should be called first

            var clusters        = new List<Ms1FeatureCluster>();
            var tempEnvelope = new double[TheoreticalEnvelope.Count];
            var tempEnvelope2 = new double[TheoreticalEnvelope.Count];
            var ms1ScanNums = Run.GetMs1ScanVector();
            foreach (var seedCell in GetSeedCells())
            {
                if (_checkedOut[seedCell.Row][seedCell.Col]) continue;
                var seedMass    = AccurateMass[seedCell.Row][seedCell.Col];
                var seedMz      = seedCell.Peaks[seedCell.RefIsotopeInternalIndex].Mz;
                var seedCharge  = ChargeRange.Min + seedCell.Row;
                var seedScanNum = ms1ScanNums[seedCell.Col];
                
                var massTol     = MzTolerance.GetToleranceAsTh(seedMass);
                var newCluster = new Ms1FeatureCluster(Run, (byte)ChargeRange.Min, TheoreticalEnvelope, seedMass, seedCharge, seedMz, seedScanNum);

                Array.Clear(tempEnvelope2, 0, tempEnvelope2.Length);
                seedCell.Peaks.SumEnvelopeTo(tempEnvelope2);
                newCluster.AddMember(seedCell);
                var normalizedElutionLength = newCluster.NetLength;

                var neighbors = new Queue<ObservedEnvelope>();

                neighbors.Enqueue(seedCell); // pick a seed
                _checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    var charge = cell.Row + ChargeRange.Min;
                    var chargeNeighborGap = GetScatteredChargeLength(charge);
                    var minRw = Math.Max(cell.Row - chargeNeighborGap, Rows.First());
                    var maxRw = Math.Min(cell.Row + chargeNeighborGap, Rows.Last());
                 
                    for (var k = 0; k < 5; k++)
                    {
                        var j = cell.Col;
                        if (k < 3) j += k;
                        else j -= (k - 2);
                        
                        if (j < Cols.First() || j > Cols.Last()) continue;

                        for (var i = minRw; i <= maxRw; i++)
                        {
                            if (_checkedOut[i][j]) continue;
                            if (!(AccurateMass[i][j] > 0)) continue;
                            if (Math.Abs(seedMass - AccurateMass[i][j]) > massTol) continue;

                            if (DistanceMap[i][j] > GetBcDistTh(normalizedElutionLength) || CorrelationMap[i][j] < GetCorrTh(normalizedElutionLength)) continue;
                            
                            // Summing envelope from paired isotopic envelopes 
                            Array.Clear(tempEnvelope, 0, TheoreticalEnvelope.Count);
                            cell.Peaks.SumEnvelopeTo(tempEnvelope);
                            FeatureMatrix[i][j].SumEnvelopeTo(tempEnvelope);

                            var newDivergence = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                            if ((DistanceMap[i][j] < 0.04 && CorrelationMap[i][j] > 0.9)|| (newDivergence < DistanceMap[i][j] && newDivergence < DistanceMap[cell.Row][cell.Col]))
                            {
                                var newEnvelope = new ObservedEnvelope(i, j, FeatureMatrix[i][j], TheoreticalEnvelope);
                                neighbors.Enqueue(newEnvelope);
                                newCluster.AddMember(newEnvelope);
                                _checkedOut[i][j] = true;
                                FeatureMatrix[i][j].SumEnvelopeTo(tempEnvelope2);
                                normalizedElutionLength = newCluster.NetLength;
                            }
                        }
                    }
                }

                for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) _checkedOut[i][j] = true;

                var summedCorr = TheoreticalEnvelope.GetPearsonCorrelation(tempEnvelope2);
                var summedBc = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope2);
                if (summedCorr < 0.2 || summedBc > 0.2) continue;

                EvaluateCluster(ref newCluster);

                if (newCluster.Envelopes.Count < 1) continue;
                if (newCluster.GoodEnvelopeCount < 1) continue;
                if (filtering && !newCluster.GoodEnough) continue;

                clusters.Add(newCluster);
            }

            return clusters;
        }

        protected virtual void EvaluateCluster(ref Ms1FeatureCluster cluster)
        {
            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;
            var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);

            cluster.ClearMember();
            for (var col = minCol; col <= maxCol; col++)
            {
                for (var row = minRow; row <= maxRow; row++)
                {
                    var mass = AccurateMass[row][col];

                    if (mass > 0 && Math.Abs(cluster.RepresentativeMass - mass) < massTol)
                    {
                        var obsEnv = new ObservedEnvelope(row, col, FeatureMatrix[row][col], TheoreticalEnvelope);
                        cluster.AddMember(obsEnv);
                    }
                    else 
                    {
                        var observedPeaks = Spectrums[col].GetAllIsotopePeaks(cluster.RepresentativeMass, row + ChargeRange.Min, TheoreticalEnvelope, MzTolerance);
                        var obsEnv = new ObservedEnvelope(row, col, observedPeaks, TheoreticalEnvelope);
                        if (obsEnv.NumberOfPeaks < 3) continue;
                        cluster.AddMember(obsEnv);    
                    }
                }
            }

            if (cluster.Envelopes.Count < 1) return;
            
            cluster.UpdateScores(Spectrums);
            if (cluster.GoodEnvelopeCount > 0)
            {
                CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
                cluster.ExpandElutionRange();
            }
        }

        private const int MinXicWindowLength = 10;
        protected void CalculateXicCorrelationOverTimeBetweenIsotopes(Ms1FeatureCluster cluster)
        {
            var maxCol = cluster.MaxCol;
            var minCol = cluster.MinCol;
            var maxRow = cluster.MaxRow;
            var minRow = cluster.MinRow;
            var colLen = maxCol - minCol + 1;

            if (colLen < MinXicWindowLength)
            {
                minCol = Math.Max(minCol - (int)((MinXicWindowLength - colLen) * 0.5), 0);

                if (minCol == 0)
                {
                    maxCol = Math.Min(minCol + MinXicWindowLength - 1, NScans - 1);
                }
                else
                {
                    maxCol = Math.Min(maxCol + (int)((MinXicWindowLength - colLen) * 0.5), NScans - 1);
                    if (maxCol == NScans - 1) minCol = Math.Max(maxCol - MinXicWindowLength + 1, 0);
                }
                colLen = maxCol - minCol + 1;
            }

            var xicProfile = new double[TheoreticalEnvelope.NhighAbundantIsotopes][];
            for (var i = 0; i < xicProfile.Length; i++) xicProfile[i] = new double[colLen];

            var monoMass = cluster.RepresentativeMass;
            var massTol = MzTolerance.GetToleranceAsTh(monoMass);
            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    var mass = AccurateMass[row][col];
                    if (!(mass > 0)) continue;
                    if (Math.Abs(monoMass - mass) > massTol) continue;
                    for (var k = 0; k < xicProfile.Length; k++)
                    {
                        var isoIndex = TheoreticalEnvelope.SortedIndexByIntensity[k];
                        if (FeatureMatrix[row][col][isoIndex] == null) continue;
                        xicProfile[k][col - minCol] += FeatureMatrix[row][col][isoIndex].Intensity;
                    }
                }
            }

            // smoothing
            //for (var k = 0; k < n; k++) xicProfile[k] = Smoother.Smooth(xicProfile[k]);
            /*var bestXicCorr = 0d;
            for (var i = 0; i < n; i++)
            {
                for (var j = i + 1; j < n; j++)
                {
                    var corr = FitScoreCalculator.GetPearsonCorrelation(xicProfile[i], xicProfile[j]);
                    if (corr > bestXicCorr) bestXicCorr = corr;
                }
            }
            return bestXicCorr;*/
            var temp = new List<double>();
            for (var i = 1; i < xicProfile.Length; i++)
            {
                var corr = FitScoreCalculator.GetBhattacharyyaDistance(xicProfile[0], xicProfile[i]);
                if (double.IsNaN(corr)) continue;
                temp.Add(corr);
                //if (corr > bestXicCorr) bestXicCorr = corr;
            }

            //cluster.ClusteringScore = temp.Median();
            //cluster.ClusteringScore2 = temp.Min();
            if (temp.Count > 0)
            {
                cluster.SetScore(Ms1FeatureScore.XicCorrMean, temp.Median());
                cluster.SetScore(Ms1FeatureScore.XicCorrMin, temp.Min());
            }
            else
            {
                cluster.SetScore(Ms1FeatureScore.XicCorrMean, 1);
                cluster.SetScore(Ms1FeatureScore.XicCorrMin, 1);
            }
            //return temp.Mean();
        }

        internal static double CalculatePoissonScore(int n, int k, int n1, int k1)
        {
            var lambda = ((double)n1 / (double)n) * k;
            var pvalue = 1 - Poisson.CDF(lambda, k1);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }
    }

    public class IsotopeList : List<Isotope>
    {
        public double MonoIsotopeMass { get; private set; }
        public readonly double[] Envelope;
        public readonly double[] EnvelopePdf;
        public readonly byte[] SortedIndexByIntensity;
        public readonly byte[] IsotopeRanking;
        public readonly int NhighAbundantIsotopes;
        
        public IsotopeList(double mass, int maxNumOfIsotopes, double relativeIntensityThreshold = 0.1)
        {
            MonoIsotopeMass = mass;
            var isoEnv = Averagine.GetIsotopomerEnvelope(MonoIsotopeMass);
            var isotopeRankings = ArrayUtil.GetRankings(isoEnv.Envolope);

            for (var i = 0; i < isoEnv.Envolope.Length; i++)
            {
                if (isoEnv.Envolope[i] < relativeIntensityThreshold || isotopeRankings[i] > maxNumOfIsotopes) continue;

                Add(new Isotope(i, isoEnv.Envolope[i]));
            }

            Envelope = this.Select(iso => iso.Ratio).ToArray();
            SortedIndexByIntensity = new byte[Count];
            EnvelopePdf = new double[Count];
            IsotopeRanking = new byte[Count];
            var s = Envelope.Sum();
            NhighAbundantIsotopes = 0;
            for (var i = 0; i < Count; i++)
            {
                IsotopeRanking[i] = (byte)isotopeRankings[this[i].Index];
                var rankingIndex = isotopeRankings[this[i].Index] - 1;
                SortedIndexByIntensity[rankingIndex] = (byte)i;
                EnvelopePdf[i] = Envelope[i] / s;

                if (Envelope[i] > Ms1Spectrum.RelativeSignificantIntesnityThreshold) NhighAbundantIsotopes++;
            }

            NhighAbundantIsotopes = Math.Min(Math.Max(NhighAbundantIsotopes, 4), Count);
        }

        public Isotope GetIsotopeRankedAt(int ranking)
        {
            return this[SortedIndexByIntensity[ranking - 1]];
        }

        public Isotope GetMostAbundantIsotope()
        {
            return GetIsotopeRankedAt(1);
        }
        
        public double GetPearsonCorrelation(double[] observedIsotopeEnvelop)
        {
            var m1 = 0.0;
            var m2 = 0.0;

            for (var i = 0; i < Count; i++)
            {
                m1 += Envelope[i];
                if (observedIsotopeEnvelop[i] > 0)
                {
                    m2 += observedIsotopeEnvelop[i];
                }
            }

            m1 /= Count;
            m2 /= Count;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < Count; i++)
            {
                var d1 = Envelope[i] - m1;
                var d2 = observedIsotopeEnvelop[i] - m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0f : cov / Math.Sqrt(s1 * s2);
        }
        
        public double GetBhattacharyyaDistance(double[] observedIsotopeEnvelop)
        {
            return FitScoreCalculator.GetBhattacharyyaDistance(Envelope, observedIsotopeEnvelop, Envelope.Length);
        }
    }

}