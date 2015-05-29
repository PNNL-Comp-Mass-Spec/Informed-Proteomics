using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Principal;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakMatrix
    {
        public LcMsPeakMatrix(LcMsRun run, int minScanCharge, int maxScanCharge, int maxThreadCount)
        {
            Run = run;
            MinScanCharge = minScanCharge;
            MaxScanCharge = maxScanCharge;
            _maxThreadCount = maxThreadCount;

            _ms1PeakList = new List<Ms1Peak>();
            Ms1Spectra = new List<Ms1Spectrum>();
            var ms1ScanNums = run.GetMs1ScanVector();

            NColumns = ms1ScanNums.Length;
            NRows = maxScanCharge - minScanCharge + 1;
            for (var i = 0; i < Math.Min(ms1ScanNums.Length, ushort.MaxValue); i++)
            {
                var ms1Spec = run.GetMs1Spectrum(ms1ScanNums[i]);
                Ms1Spectra.Add(ms1Spec);
                _ms1PeakList.AddRange((Ms1Peak[])ms1Spec.Peaks);
            }
            _ms1PeakList.Sort();
            
            MzTolerance = new Tolerance(5);
            _featureMatrixCreated = false;

            _distProfileAcrossTime = new double[NColumns];
            _distProfileAcrossCharge = new double[NRows];
            _xic = new double[NColumns];

            _comparer = new MzComparerWithBinning(27);
        }

        public IList<LcMsPeakCluster> FindFeatures(double targetMass)
        {
            SetTargetMass(targetMass);
            var features = ClusterMs1Peaks(true);
            return features;
        }

        public Ms1Spectrum GetSpectrum(int ms1ScanNum)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var idx = Array.BinarySearch(ms1ScanNums, ms1ScanNum);
            return idx < 0 ? null : Ms1Spectra[idx];
        }

        public LcMsRun Run;
        public readonly List<Ms1Spectrum> Ms1Spectra;
        public readonly int MinScanCharge;
        public readonly int MaxScanCharge;
        public readonly Tolerance MzTolerance;

        protected readonly int NColumns;
        protected readonly int NRows;

        private readonly int _maxThreadCount;
        private readonly MzComparerWithBinning _comparer;

        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;

        private readonly double[] _distProfileAcrossCharge;
        private readonly double[] _distProfileAcrossTime;
        private readonly double[] _xic;

        private bool _featureMatrixCreated;
        protected bool[][] CheckedOut;
        protected double[][] CorrelationMap;
        protected double[][] DistanceMap;
        protected double[][] AccurateMass;
        protected Ms1Peak[][][] FeatureMatrix;

        protected TheoreticalIsotopeEnvelope TheoreticalEnvelope;
        protected int[] Rows;
        protected int[] Cols;

        protected double TargetMass;
        protected bool IsAccurateMass;
        
        private void InitFeatureMatrix()
        {
            CorrelationMap = new double[NRows][];
            DistanceMap = new double[NRows][];
            AccurateMass = new double[NRows][];
            FeatureMatrix = new Ms1Peak[NRows][][];
            CheckedOut = new bool[NRows][];

            for (var i = 0; i < NRows; i++)
            {
                CheckedOut[i] = new bool[NColumns];
                CorrelationMap[i] = new double[NColumns];
                FeatureMatrix[i] = new Ms1Peak[NColumns][];
                DistanceMap[i] = new double[NColumns];
                AccurateMass[i] = new double[NColumns];

                for (var j = 0; j < NColumns; j++)
                {
                    FeatureMatrix[i][j] = new Ms1Peak[MaxEnvelopeLength];
                }
            }
            _featureMatrixCreated = true;
        }
        
        protected void SetTargetMass(double targetMass, bool isAccurateMass = false)
        {
            TargetMass = targetMass;
            TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, 0.1d);

            var chargeLowerBound = Math.Floor(TargetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Ceiling(TargetMass / Run.MinMs1Mz);

            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);

            Rows = Enumerable.Range((int) rowLb, (int) (rowUb - rowLb + 1)).ToArray();

            IsAccurateMass = isAccurateMass;
        }
       
        protected void BuildFeatureMatrix()
        {
            if (!_featureMatrixCreated) InitFeatureMatrix();

            var queryMassBinNum = _comparer.GetBinNumber(TargetMass);

            var options             = new ParallelOptions();
            var observedRows        = new bool[NRows];
            var observedCols        = new bool[NColumns];

            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            var mostAbuIsotopeInternalIndex = TheoreticalEnvelope.IndexOrderByRanking[0];
            var mostAbuIsotopeIndex = TheoreticalEnvelope.GetMostAbundantIsotope().Index;
            
            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;

            var signalToNoiseRatioCutoff = (IsAccurateMass) ? 0 : 1.4826;
            var nPeaksCutoff = (TargetMass > 2000) ? 3 : 2;

            Parallel.ForEach(Rows, options, row =>
            {
                Array.Clear(CorrelationMap[row], 0, NColumns);
                Array.Clear(CheckedOut[row], 0, NColumns);
                Array.Clear(DistanceMap[row], 0, NColumns);
                Array.Clear(AccurateMass[row], 0, NColumns);

                for (var col = 0; col < NColumns; col++)
                {
                    Array.Clear(FeatureMatrix[row][col], 0, FeatureMatrix[row][col].Length);
                }

                var charge = row + MinScanCharge;
                for (var k = 0; k < TheoreticalEnvelope.Size; k++)
                {
                    var i = TheoreticalEnvelope.IndexOrderByRanking[k];
                    var isotopeIndex = TheoreticalEnvelope.Isotopes[i].Index;

                    double isotopeMzLb;
                    double isotopeMzUb;

                    if (IsAccurateMass)
                    {
                        var isotopeMz = Ion.GetIsotopeMz(TargetMass, charge, isotopeIndex);
                        var mzTol = MzTolerance.GetToleranceAsTh(isotopeMz);
                        isotopeMzLb = isotopeMz - mzTol;
                        isotopeMzUb = isotopeMz + mzTol;
                    }
                    else
                    {
                        isotopeMzLb = (k == 0) ? Ion.GetIsotopeMz(_comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(_comparer.GetMzAverage(queryMassBinNum - 1), charge, isotopeIndex);
                        isotopeMzUb = (k == 0) ? Ion.GetIsotopeMz(_comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(_comparer.GetMzAverage(queryMassBinNum + 1), charge, isotopeIndex);
                    }

                    if (isotopeMzLb < minMs1Mz || isotopeMzUb > maxMs1Mz) continue;
                    var st = _ms1PeakList.BinarySearch(new Ms1Peak(isotopeMzLb, 0, 0));

                    if (st < 0) st = ~st;

                    for (var j = st; j < _ms1PeakList.Count; j++)
                    {
                        var ms1Peak         = _ms1PeakList[j];
                        if (ms1Peak.Mz > isotopeMzUb) break;
                        var col = ms1Peak.Ms1SpecIndex;

                        if (k == 0) // most abundant peak
                        {
                            if (FeatureMatrix[row][col][i] == null || ms1Peak.Intensity > FeatureMatrix[row][col][i].Intensity)
                            {
                                FeatureMatrix[row][col][i] = ms1Peak;
                                AccurateMass[row][col] = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                            }
                        }
                        else
                        {
                            if (IsAccurateMass)
                            {
                                if (!(AccurateMass[row][col] > 0)) AccurateMass[row][col] = TargetMass;
                            }
                            else
                            {
                                // mass binning search mode, most abundant peak should be exist
                                var mostAbuPeak = FeatureMatrix[row][col][mostAbuIsotopeInternalIndex];
                                if (mostAbuPeak == null) continue;
                                var expectedPeakMz = mostAbuPeak.Mz +
                                                     (Constants.C13MinusC12 * (isotopeIndex - mostAbuIsotopeIndex)) / charge;
                                if (Math.Abs(expectedPeakMz - ms1Peak.Mz) > MzTolerance.GetToleranceAsTh(ms1Peak.Mz))
                                    continue;                                
                            }

                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (FeatureMatrix[row][col][i] != null)
                            {
                                var tmpPeak = FeatureMatrix[row][col][i];
                                var bc1 = TheoreticalEnvelope.GetBhattacharyyaDistance(FeatureMatrix[row][col]);
                                FeatureMatrix[row][col][i] = ms1Peak;
                                var bc2 = TheoreticalEnvelope.GetBhattacharyyaDistance(FeatureMatrix[row][col]);
                                if (bc1 < bc2) FeatureMatrix[row][col][i] = tmpPeak;
                            }
                            else
                            {
                                FeatureMatrix[row][col][i] = ms1Peak;

                            }
                        }
                    }
                }

                for (var col = 0; col < NColumns; col++)
                {
                    if (!(AccurateMass[row][col] > 0)) continue;

                    var nPeaks = FeatureMatrix[row][col].Count(p => p != null && p.Active);
                    if (IsAccurateMass)
                    {
                        if (nPeaks >= nPeaksCutoff)
                        {
                            CorrelationMap[row][col] = TheoreticalEnvelope.GetPearsonCorrelation(FeatureMatrix[row][col]);
                            DistanceMap[row][col] = TheoreticalEnvelope.GetBhattacharyyaDistance(FeatureMatrix[row][col]);
                            if (!observedRows[row]) observedRows[row] = true;
                            if (!observedCols[col]) observedCols[col] = true;
                        }
                        else
                        {
                            AccurateMass[row][col] = 0d;
                        }
                    }
                    else
                    {
                        var mostAbuPeakIntensity = FeatureMatrix[row][col][mostAbuIsotopeInternalIndex].Intensity;
                        var signalToNoiseRatio = mostAbuPeakIntensity / Ms1Spectra[col].MedianIntensity;
                        if (signalToNoiseRatio > signalToNoiseRatioCutoff && nPeaks >= nPeaksCutoff)
                        {
                            CorrelationMap[row][col] = TheoreticalEnvelope.GetPearsonCorrelation(FeatureMatrix[row][col]);
                            DistanceMap[row][col] = TheoreticalEnvelope.GetBhattacharyyaDistance(FeatureMatrix[row][col]);
                            if (!observedRows[row]) observedRows[row] = true;
                            if (!observedCols[col]) observedCols[col] = true;
                        }
                        else
                        {
                            AccurateMass[row][col] = 0d;
                        }
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

        private IEnumerable<ObservedIsotopeEnvelope> GetSeedCells()
        {
            var seedList = new List<KeyValuePair<double, ObservedIsotopeEnvelope>>();
            var mostAbuIsotopeInternalIndex = TheoreticalEnvelope.IndexOrderByRanking[0];
            var ms1ScanNums = Run.GetMs1ScanVector();

            var nPeaksCutoff = 2; 
            if (TargetMass > 25000) nPeaksCutoff = 5;
            else if (TargetMass > 8000) nPeaksCutoff = 4;
            else if (TargetMass > 2000) nPeaksCutoff = 3;
            else nPeaksCutoff = 2;

            foreach (var i in Rows)
            {
                foreach (var j in Cols)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;

                    var mostAbuPeakIntensity = FeatureMatrix[i][j][mostAbuIsotopeInternalIndex].Intensity;
                    var signalToNoiseRatio = mostAbuPeakIntensity / Ms1Spectra[j].MedianIntensity;
                    if (signalToNoiseRatio < 3) continue;

                    if (FeatureMatrix[i][j].Count(p => p != null && p.Active) < nPeaksCutoff) continue;

                    var bcDist = DistanceMap[i][j];
                    var corr = CorrelationMap[i][j];

                    if (bcDist > 0.3 && corr < 0.5) continue;

                    var seed = new ObservedIsotopeEnvelope(TargetMass, i + MinScanCharge, ms1ScanNums[j], FeatureMatrix[i][j], TheoreticalEnvelope);
                    seedList.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(bcDist, seed));
                }
            }
            return seedList.OrderBy(x => x.Key).Select(x => x.Value);
        }
        
        private List<LcMsPeakCluster> ClusterMs1Peaks(bool filtering = false)
        {
            BuildFeatureMatrix(); // should be called first

            var clusters        = new List<LcMsPeakCluster>();
            var tempEnvelope = new double[TheoreticalEnvelope.Size];
            var tempEnvelope2 = new double[TheoreticalEnvelope.Size];
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            foreach (var seed in GetSeedCells())
            {
                var row = seed.Charge - MinScanCharge;
                var col = ms1ScanNumToIndex[seed.ScanNum];

                if (CheckedOut[row][col]) continue;
                var seedMass = AccurateMass[row][col];
                var massTol = MzTolerance.GetToleranceAsTh(seedMass);
                var newCluster = new LcMsPeakCluster(Run, seed);

                Array.Clear(tempEnvelope2, 0, tempEnvelope2.Length);
                seed.Peaks.SumEnvelopeTo(tempEnvelope2);

                var neighbors = new Queue<ObservedIsotopeEnvelope>();
                neighbors.Enqueue(seed); // pick a seed
                CheckedOut[row][col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    var charge = cell.Charge;
                    var chargeNeighborGap = GetScatteredChargeLength(charge);
                    var minRw = Math.Max(charge - MinScanCharge - chargeNeighborGap, Rows.First());
                    var maxRw = Math.Min(charge - MinScanCharge + chargeNeighborGap, Rows.Last());
                    
                    var currRow = charge - MinScanCharge;
                    var currCol = ms1ScanNumToIndex[cell.ScanNum];
                 
                    for (var k = 0; k < 5; k++)
                    {
                        var j = currCol;
                        if (k < 3) j += k;
                        else j -= (k - 2);
                        
                        if (j < Cols.First() || j > Cols.Last()) continue;

                        for (var i = minRw; i <= maxRw; i++)
                        {
                            if (CheckedOut[i][j]) continue;
                            if (!(AccurateMass[i][j] > 0)) continue;
                            if (Math.Abs(seedMass - AccurateMass[i][j]) > massTol) continue;

                            if (DistanceMap[i][j] > GetBcDistTh(newCluster.NetLength) || CorrelationMap[i][j] < GetCorrTh(newCluster.NetLength)) continue;
                            
                            // Summing envelope from paired isotopic envelopes 
                            Array.Clear(tempEnvelope, 0, tempEnvelope.Length);
                            cell.Peaks.SumEnvelopeTo(tempEnvelope);
                            FeatureMatrix[i][j].SumEnvelopeTo(tempEnvelope);

                            var newDivergence = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                            if ((DistanceMap[i][j] < 0.04 && CorrelationMap[i][j] > 0.9)|| (newDivergence < DistanceMap[i][j] && newDivergence < DistanceMap[currRow][currCol]))
                            {
                                var envelope = new ObservedIsotopeEnvelope(AccurateMass[i][j], i + MinScanCharge, ms1ScanNums[j], FeatureMatrix[i][j], TheoreticalEnvelope);

                                neighbors.Enqueue(envelope);
                                newCluster.AddObservedEnvelope(envelope);
                                CheckedOut[i][j] = true;
                                FeatureMatrix[i][j].SumEnvelopeTo(tempEnvelope2);
                            }
                        }
                    }
                }

                for (var i = newCluster.MinCharge - MinScanCharge; i <= newCluster.MaxCharge - MinScanCharge; i++)
                    for (var j = ms1ScanNumToIndex[newCluster.MinScanNum]; j <= ms1ScanNumToIndex[newCluster.MaxScanNum]; j++) CheckedOut[i][j] = true;

                //var summedCorr = TheoreticalEnvelope.GetPearsonCorrelation(tempEnvelope2);
                //var summedBc = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope2);
                //if (summedCorr < 0.2 || summedBc > 0.2) continue;

                /////////////////////////////////////////////// initial evaluation ////////////////////////////////////
                //EvaluateCluster(ref newCluster);

                var cluster = FindLcMsPeakCluster(newCluster.RepresentativeMass, newCluster.RepresentativeScanNum, newCluster.RepresentativeCharge, newCluster);

                //if (newCluster.Envelopes.Count < 1) continue;
                //if (newCluster.GoodEnvelopeCount < 1) continue;
                //if (filtering && !newCluster.GoodEnough) continue;
                /////////////////////////////////////////////// initial evaluation ////////////////////////////////////

                clusters.Add(newCluster);
            }
            return clusters;
        }



        private double GetBcDistTh(double normalizedElutionLen)
        {
            if (TargetMass < 15000)
            {
                if (normalizedElutionLen < 0.005) return 0.6;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.2;
                return 0.1;
            }
            else if (TargetMass < 25000)
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

        private double GetCorrTh(double normalizedElutionLen)
        {
            if (TargetMass < 15000)
            {
                if (normalizedElutionLen < 0.005) return 0.3;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.6;
                return 0.7;
            }
            else if (TargetMass < 25000)
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

        private int GetScatteredChargeLength(int charge)
        {
            var scatteredLength = 2;
            if (charge > 40) scatteredLength = 12;
            else if (charge > 30) scatteredLength = 9;
            else if (charge > 20) scatteredLength = 6;
            else if (charge > 10) scatteredLength = 3;
            return scatteredLength;
        }

        private Tuple<int, int> GetElutionWindow(int col, double halfWindowSize, int minScanCount = 5)
        {
            var minCol = col;
            var maxCol = col;
            var ms1ScanNums = Run.GetMs1ScanVector();

            var elutionTime = Run.GetElutionTime(ms1ScanNums[col]);

            for (var j = col - 1; j >= Cols.First(); j--)
            {
                //if (j < Cols.First() || j > Cols.Last()) break;
                if (elutionTime - Run.GetElutionTime(ms1ScanNums[j]) > halfWindowSize && (col - j) >= minScanCount) break;
                if (j < minCol) minCol = j;
            }

            for (var j = col + 1; j < Cols.Last(); j++)
            {
                //if (j < Cols.First() || j > Cols.Last()) continue;
                if (Run.GetElutionTime(ms1ScanNums[j]) - elutionTime > halfWindowSize && (j - col) >= minScanCount) break;
                if (j > maxCol) maxCol = j;
            }
            return new Tuple<int, int>(minCol, maxCol);
        }
        
        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);

        public LcMsPeakCluster FindLcMsPeakCluster_old(double targetMass, int targetScanNum, int targetCharge, LcMsFeature initialFeature = null)
        {
            SetTargetMass(targetMass, true);

            BuildFeatureMatrix(); // should be called first

            const double initialElutionPeriod = 0.5; // 1 minute
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            var minCol = 0;
            var maxCol = 0;
            var row = targetCharge - MinScanCharge;
            var col = Array.BinarySearch(ms1ScanNums, targetScanNum);
            if (col < 0) col = ~col;

            var bcCutoff = OtsuThreshold.GetThreshold(DistanceMap, 0, 0.5, 0.005, Rows.First(), Rows.Last(), Cols.First(), Cols.Last(), 1e-9, 0.5);

            //Console.WriteLine("-------------------------"); Console.WriteLine("cutoff = {0}", bcCutoff);

            // get seed 
            var tempWindow = GetElutionWindow(col, initialElutionPeriod);
            minCol = tempWindow.Item1;
            maxCol = tempWindow.Item2;
            //var minRow = Math.Max(row - 5, Rows.First());
            //var maxRow = Math.Min(row + 5, Rows.Last());

            var seeds = new List<KeyValuePair<double, Tuple<int, int>>>();
            for (var j = minCol; j <= maxCol; j++)
            {
                //for (var i = minRow; i <= maxRow; i++)
                foreach(var i in Rows)
                {
                    if (AccurateMass[i][j] > 0 && DistanceMap[i][j] < bcCutoff)
                    {
                        var dist = Math.Abs(col - j);
                        seeds.Add(new KeyValuePair<double, Tuple<int, int>>(dist, new Tuple<int, int>(i, j)));
                    }
                }
            }



            LcMsPeakCluster feature = null;
            foreach (var pt in seeds.OrderBy(x => x.Key).Select(x => x.Value))
            {
                var i = pt.Item1;
                var j = pt.Item2;

                //Console.WriteLine("seed = {0}, {1}", i+MinScanCharge, ms1ScanNums[j]);

                feature = ClusterMs1Peaks(i, j, bcCutoff);

                break;
            }

            if (feature == null)
            {
                return new LcMsPeakCluster(Run, TargetMass, targetCharge, 0, ms1ScanNums[col], 0d);
            }

            /*
            if (initialFeature != null)
            {
                minCol = ms1ScanNumToIndex[initialFeature.MinScanNum];
                maxCol = ms1ScanNumToIndex[initialFeature.MaxScanNum];
            }
            else
            {
                var tempWindow = GetElutionWindow(col, initialElutionPeriod);
                minCol = tempWindow.Item1;
                maxCol = tempWindow.Item2;
            }

            // 1) determine charge state range
            var minMaxRow = FindMinMaxRow(row, minCol, maxCol);
            var minRow = minMaxRow.Item1;
            var maxRow = minMaxRow.Item2;

            // 2) determine elution period
            var minMaxCol = FindMinMaxCol(col, minRow, maxRow, minCol, maxCol);
            minCol = minMaxCol.Item1;
            maxCol = minMaxCol.Item2;

            // 3) collect envelopes
            var feature = GetLcMsPeakCluster(minRow, maxRow, minCol, maxCol);

            if (feature == null)
            {
                return new LcMsPeakCluster(Run, TargetMass, targetCharge, 0, ms1ScanNums[col], 0d);
            }

            feature.EnvelopeDistanceScoreAcrossCharge = new double[maxRow - minRow + 1];
            feature.EnvelopeDistanceScoreAcrossTime = new double[maxCol - minCol + 1];

            Array.Copy(_distProfileAcrossTime, minCol, feature.EnvelopeDistanceScoreAcrossTime, 0, feature.EnvelopeDistanceScoreAcrossTime.Length);
            Array.Copy(_distProfileAcrossCharge, minRow, feature.EnvelopeDistanceScoreAcrossCharge, 0, feature.EnvelopeDistanceScoreAcrossCharge.Length);
            */
            

            // 4) determine abundance
            Array.Clear(_xic, 0, _xic.Length);
            foreach (var envelope in feature.Envelopes)
            {
                var envCol = ms1ScanNumToIndex[envelope.ScanNum];
                _xic[envCol] += envelope.Abundance;
            }

            var smoothedXic = Smoother.Smooth(_xic);
            var abundance = 0d;
            for (var k = 0; k < smoothedXic.Length - 1; k++)
            {
                var centerIntensity = 0.5 * (Math.Max(0, smoothedXic[k]) + Math.Max(0, smoothedXic[k + 1]));

                if (!(centerIntensity > 0)) continue;
                var timeInterval = Run.GetElutionTime(ms1ScanNums[k + 1]) - Run.GetElutionTime(ms1ScanNums[k]);

                abundance += centerIntensity * timeInterval;
            }
            feature.SetAbundance(abundance);

            return feature;
        }

        private LcMsPeakCluster ClusterMs1Peaks(int row, int col, double bcDistCutoff)
        {
            if (CheckedOut[row][col]) return null;
            var seedMass = AccurateMass[row][col];
            //var massTol = MzTolerance.GetToleranceAsTh(seedMass);
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var seed = new ObservedIsotopeEnvelope(seedMass, row + MinScanCharge, ms1ScanNums[col], FeatureMatrix[row][col], TheoreticalEnvelope);
            
            var newCluster = new LcMsPeakCluster(Run, seed);
            var neighbors = new Queue<ObservedIsotopeEnvelope>();
            neighbors.Enqueue(seed); // pick a seed
            CheckedOut[row][col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    var charge = cell.Charge;
                    var chargeNeighborGap = GetScatteredChargeLength(charge);
                    var minRw = Math.Max(charge - MinScanCharge - chargeNeighborGap, Rows.First());
                    var maxRw = Math.Min(charge - MinScanCharge + chargeNeighborGap, Rows.Last());

                    //var currRow = charge - MinScanCharge;
                    var currCol = ms1ScanNumToIndex[cell.ScanNum];

                    for (var k = 0; k < 5; k++)
                    {
                        var j = currCol;
                        if (k < 3) j += k;
                        else j -= (k - 2);

                        if (j < Cols.First() || j > Cols.Last()) continue;

                        for (var i = minRw; i <= maxRw; i++)
                        {
                            if (CheckedOut[i][j]) continue;
                            if (!(AccurateMass[i][j] > 0)) continue;
                            //if (Math.Abs(seedMass - AccurateMass[i][j]) > massTol) continue;
                            //if (DistanceMap[i][j] > GetBcDistTh(newCluster.NetLength) || CorrelationMap[i][j] < GetCorrTh(newCluster.NetLength)) continue;
                            //if ((DistanceMap[i][j] < 0.04 && CorrelationMap[i][j] > 0.9) || (newDivergence < DistanceMap[i][j] && newDivergence < DistanceMap[currRow][currCol]))
                            if (DistanceMap[i][j] < bcDistCutoff)
                            {
                                var envelope = new ObservedIsotopeEnvelope(AccurateMass[i][j], i + MinScanCharge, ms1ScanNums[j], FeatureMatrix[i][j], TheoreticalEnvelope);
                                neighbors.Enqueue(envelope);
                                newCluster.AddObservedEnvelope(envelope);
                                CheckedOut[i][j] = true;
                                //FeatureMatrix[i][j].SumEnvelopeTo(tempEnvelope2);
                            }
                        }
                    }
                }

                for (var i = newCluster.MinCharge - MinScanCharge; i <= newCluster.MaxCharge - MinScanCharge; i++)
                    for (var j = ms1ScanNumToIndex[newCluster.MinScanNum]; j <= ms1ScanNumToIndex[newCluster.MaxScanNum]; j++) CheckedOut[i][j] = true;

                //var cluster = FindLcMsPeakCluster(newCluster.RepresentativeMass, newCluster.RepresentativeScanNum, newCluster.RepresentativeCharge, newCluster);

                //if (newCluster.Envelopes.Count < 1) continue;
                //if (newCluster.GoodEnvelopeCount < 1) continue;
                //if (filtering && !newCluster.GoodEnough) continue;
                /////////////////////////////////////////////// initial evaluation ////////////////////////////////////

                //clusters.Add(newCluster);
            //}
            //return clusters;
            return newCluster;
        }


        public LcMsPeakCluster FindLcMsPeakCluster(double targetMass, int targetScanNum, int targetCharge, LcMsFeature initialFeature = null)
        {
            SetTargetMass(targetMass, true);

            BuildFeatureMatrix(); // should be called first

            const double initialElutionPeriod = 0.5; // 1 minute
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            var minCol = 0;
            var maxCol = 0;
            var row = targetCharge - MinScanCharge;
            var col = Array.BinarySearch(ms1ScanNums, targetScanNum);
            if (col < 0) col = ~col;
            
            if (initialFeature != null)
            {
                minCol = ms1ScanNumToIndex[initialFeature.MinScanNum];
                maxCol = ms1ScanNumToIndex[initialFeature.MaxScanNum];
            }
            else
            {
                var tempWindow = GetElutionWindow(col, initialElutionPeriod);
                minCol = tempWindow.Item1;
                maxCol = tempWindow.Item2;

                //Console.WriteLine("initial search scans : {0} - {1}", ms1ScanNums[minCol], ms1ScanNums[maxCol]);
            }

            // 1) determine charge state range
            var minMaxRow = FindMinMaxRow(row, minCol, maxCol);
            var minRow = minMaxRow.Item1;
            var maxRow = minMaxRow.Item2;

            //Console.WriteLine("{0} - {1}", minRow + MinScanCharge, maxRow + MinScanCharge);

            // 2) determine elution period
            var minMaxCol = FindMinMaxCol(col, minRow, maxRow, minCol, maxCol);
            minCol = minMaxCol.Item1;
            maxCol = minMaxCol.Item2;

            //Console.WriteLine("{0} - {1}", ms1ScanNums[minCol], ms1ScanNums[maxCol]);
            
            // 3) collect envelopes
            var feature = GetLcMsPeakCluster(minRow, maxRow, minCol, maxCol);

            if (feature == null)
            {
                return new LcMsPeakCluster(Run, TargetMass, targetCharge, 0, ms1ScanNums[col], 0d);
            }

            feature.EnvelopeDistanceScoreAcrossCharge = new double[maxRow - minRow + 1];
            feature.EnvelopeDistanceScoreAcrossTime = new double[maxCol - minCol + 1];

            Array.Copy(_distProfileAcrossTime, minCol, feature.EnvelopeDistanceScoreAcrossTime, 0, feature.EnvelopeDistanceScoreAcrossTime.Length);
            Array.Copy(_distProfileAcrossCharge, minRow, feature.EnvelopeDistanceScoreAcrossCharge, 0, feature.EnvelopeDistanceScoreAcrossCharge.Length);
            
            // 4) determine abundance
            Array.Clear(_xic, 0, _xic.Length);
            foreach (var envelope in feature.Envelopes)
            {
                var envCol = ms1ScanNumToIndex[envelope.ScanNum];
                _xic[envCol] += envelope.Abundance;
            }

            
            for (var i = minRow; i <= maxRow; i++)
            {
                var sEnvelope = GetSummedEnvelope(i, i, minCol, maxCol);
                
                Console.Write(TheoreticalEnvelope.GetBhattacharyyaDistance(sEnvelope));
                Console.Write("\t");
                Console.Write(TheoreticalEnvelope.GetPearsonCorrelation(sEnvelope));

                for (var j = 0; j < sEnvelope.Length; j++)
                {
                    Console.Write("\t");
                    Console.Write(sEnvelope[j]);
                }
                Console.Write("\n");
            }
            


            var smoothedXic = Smoother.Smooth(_xic);
            var abundance = 0d;
            for (var k = 0; k < smoothedXic.Length - 1; k++)
            {
                var centerIntensity = 0.5 * (Math.Max(0, smoothedXic[k]) + Math.Max(0, smoothedXic[k + 1]));

                if (!(centerIntensity > 0)) continue;
                var timeInterval = Run.GetElutionTime(ms1ScanNums[k + 1]) - Run.GetElutionTime(ms1ScanNums[k]);

                abundance += centerIntensity * timeInterval;
            }
            feature.SetAbundance(abundance);
            
            return feature;
        }

        private LcMsPeakCluster GetLcMsPeakCluster(int minRow, int maxRow, int minCol, int maxCol)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();

            var envelopes = new List<ObservedIsotopeEnvelope>();
            var bcDistList = new List<double>();
            var corrList = new List<double>();
            var bestBcDist = 100d;
            ObservedIsotopeEnvelope bestEnvelope = null;

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;
                    var envelope = new ObservedIsotopeEnvelope(AccurateMass[i][j], i + MinScanCharge, ms1ScanNums[j], FeatureMatrix[i][j], TheoreticalEnvelope);
                    envelopes.Add(envelope);
                    bcDistList.Add(DistanceMap[i][j]);
                    corrList.Add(CorrelationMap[i][j]);

                    if (DistanceMap[i][j] < bestBcDist)
                    {
                        bestBcDist = DistanceMap[i][j];
                        bestEnvelope = envelope;
                    }
                }
            }

            if (bestEnvelope == null)
            {
                return null;
            }

            //return envelopes;
            
            var bcDistSample = bcDistList.Where(d => d < bcDistList.Median()).ToList();
            var meanStd = bcDistSample.MeanStandardDeviation();
            var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.05);
            bcCutoff = Math.Min(bcCutoff, 0.5);

            var corrSample = corrList.Where(d => d < corrList.Median()).ToList();
            meanStd = corrSample.MeanStandardDeviation();
            var corrCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.3);
            corrCutoff = Math.Min(bcCutoff, 0.9);

            var cluster = new LcMsPeakCluster(Run, bestEnvelope);
            for (var k = 0; k < bcDistList.Count; k++)
            {
                if (bcDistList[k] < bcCutoff || corrList[k] > corrCutoff) cluster.AddObservedEnvelope(envelopes[k]); 
            }

            return cluster;
        }


        private Tuple<int, int> FindMinMaxCol(int seedCol, int minRow, int maxRow, int colLb, int colUb)
        {
            const double searchStopBcThreshold = 0.3;
            
            var minColTemp = colLb;
            var maxColTemp = colUb;
            
            // sampling bc distances around local minimum using 1 min. window size
            //var tempWindow = GetElutionWindow(seedCol, 0.5);

            Array.Clear(_distProfileAcrossTime, 0, _distProfileAcrossTime.Length);

            for (var j = seedCol; j <= Cols.Last(); j++)
            {
                var summedBcDist = GetBestSummedEnvelope(minRow, maxRow, j, j);
                if (j > colUb && j > seedCol && summedBcDist > searchStopBcThreshold && _distProfileAcrossTime[j - 1] > searchStopBcThreshold) break;
                if (j > maxColTemp) maxColTemp = j;
                _distProfileAcrossTime[j] = summedBcDist;
            }

            for (var j = seedCol - 1; j >= Cols.First(); j--)
            {
                var summedBcDist = GetBestSummedEnvelope(minRow, maxRow, j, j);
                if (j < colLb && summedBcDist > searchStopBcThreshold && j < seedCol - 1 && _distProfileAcrossTime[j + 1] > searchStopBcThreshold) break;
                if (j < minColTemp) minColTemp = j;
                _distProfileAcrossTime[j] = summedBcDist;
            }
            
            var bcCutoff = OtsuThreshold.GetThreshold(_distProfileAcrossTime, 0, searchStopBcThreshold, 0.0025, 1e-9, searchStopBcThreshold);
            //Console.WriteLine("bc cut off = {0}", bcCutoff);

            /*bcDistSample.Sort();
            for (var i = 0; i < bcDistSample.Count; i++)
            {
                if (bcDistSample[i] > 0.01 && bcDistSample[i] < 0.1) bcDistSample2.Add(bcDistSample[i]);
            }

            var meanStd = bcDistSample2.MeanStandardDeviation();
            var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.01);
            bcCutoff = Math.Min(bcCutoff, 0.1);*/
            
            //var bcCutoff = 0.05;
            var maxCol = seedCol;
            var minCol = seedCol;
            for (var j = seedCol + 1; j <= maxColTemp; j++)
            {
                if (j == Cols.Last())
                {
                    maxCol = j;
                    break;
                }
                
                if (_distProfileAcrossTime[j] > bcCutoff && _distProfileAcrossTime[j + 1] > bcCutoff) break;
                maxCol = j;
            }
            for (var j = seedCol - 1; j >= minColTemp; j--)
            {
                if (j == Cols.First())
                {
                    minCol = j;
                    break;
                }
                if (_distProfileAcrossTime[j] > bcCutoff && _distProfileAcrossTime[j - 1] > bcCutoff) break;
                minCol = j;
            }

            return new Tuple<int, int>(minCol, maxCol);
        }

        private Tuple<int, int> FindMinMaxRow(int seedRow, int minCol, int maxCol)
        {
            const double searchStopBcThreshold = 0.3;
            var minRow = seedRow;
            var maxRow = seedRow;
            Array.Clear(_distProfileAcrossCharge, 0, _distProfileAcrossCharge.Length);

            //var bcDistSample = new List<double>();
            foreach (var i in Rows)
            {
                var summedBcDist = GetBestSummedEnvelope(i, i, minCol, maxCol);
                _distProfileAcrossCharge[i] = summedBcDist;
                //if (summedBcDist < searchStopBcThreshold) bcDistSample.Add(summedBcDist);

                //Console.WriteLine("\t{0}, {1}", i + MinScanCharge, summedBcDist);
            }
            

            //var meanStd = bcDistSample.MeanStandardDeviation();
            //var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.05);
            var bcCutoff = OtsuThreshold.GetThreshold(_distProfileAcrossCharge, 0, searchStopBcThreshold, 0.0025, 1e-9, searchStopBcThreshold);
            //Console.WriteLine("Charge bcCutoff = {0}", bcCutoff);


            foreach (var i in Rows)
            {
                if (_distProfileAcrossCharge[i] > bcCutoff) continue;
                if (i < minRow) minRow = i;
                if (i > maxRow) maxRow = i;
            }

            return new Tuple<int, int>(minRow, maxRow);
        }

        private double GetBestSummedEnvelope(int minRow, int maxRow, int minCol, int maxCol)
        {
            var summedBcDist = 10.0d;
            var cellList = new List<Tuple<double, int, int>>();
            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;
                    var bcDist = DistanceMap[i][j];
                    var cell = new Tuple<double, int, int>(bcDist, i, j);
                    cellList.Add(cell);
                }
            }

            if (cellList.Count < 1) return summedBcDist;
            if (cellList.Count == 1) return cellList[0].Item1;

            var tempEnvelope = new double[TheoreticalEnvelope.Size];
            var summedEnvelope = new double[TheoreticalEnvelope.Size];

            foreach (var cell in cellList.OrderBy(c => c.Item1))
            {
                FeatureMatrix[cell.Item2][cell.Item3].SumEnvelopeTo(tempEnvelope);
                var tempBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                if (tempBcDist > summedBcDist)
                {
                    Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                }
                else
                {
                    summedBcDist = tempBcDist;
                    Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                }
            }
            return summedBcDist;
        }

        private double[] GetSummedEnvelope(int minRow, int maxRow, int minCol, int maxCol)
        {
            var summedEnvelope = new double[TheoreticalEnvelope.Size];
            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;
                    //var bcDist = DistanceMap[i][j];
                    FeatureMatrix[i][j].SumEnvelopeTo(summedEnvelope);
                }
            }
            return summedEnvelope;
        }


        private void GetInitialEvaluation(LcMsPeakCluster cluster)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var minCol = ms1ScanNumToIndex[cluster.MinScanNum];
            var maxCol = ms1ScanNumToIndex[cluster.MaxScanNum];
            var minRow = cluster.MinCharge - MinScanCharge;
            var maxRow = cluster.MaxCharge + MaxScanCharge;

            var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);

            /*

            if (cluster.Envelopes.Count < 1) return;
            
            cluster.UpdateScores(Spectrums);
            if (cluster.GoodEnvelopeCount > 0)
            {
                CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
                cluster.ExpandElutionRange();
            }*/
        }
      
    }
}
