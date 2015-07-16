using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;


namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakMatrix
    {
        public LcMsPeakMatrix(LcMsRun run, LcMsFeatureLikelihood scorer = null, int maxThreadCount = 0)
        {
            Run = run;
            //MinScanCharge = minScanCharge;
            //MaxScanCharge = maxScanCharge;
            _maxThreadCount = maxThreadCount;

            _ms1PeakList = new List<Ms1Peak>();
            Ms1Spectra = new List<Ms1Spectrum>();
            var ms1ScanNums = run.GetMs1ScanVector();

            NColumns = ms1ScanNums.Length;
            NRows = MaxScanCharge - MinScanCharge + 1;
            for (var i = 0; i < Math.Min(ms1ScanNums.Length, ushort.MaxValue); i++)
            {
                var ms1Spec = run.GetMs1Spectrum(ms1ScanNums[i]);
                Ms1Spectra.Add(ms1Spec);
                _ms1PeakList.AddRange((Ms1Peak[])ms1Spec.Peaks);
            }
            _ms1PeakList.Sort();
            
            //MzTolerance = new Tolerance(5);

            _distProfileAcrossCharge = new double[NRows];
            _corrProfileAcrossCharge = new double[NRows];
            _summedEnvelopeColRange = new int[NRows, 2];
            
            _xic = new double[NColumns];
            _featureMatrix = null;
            Comparer = new MzComparerWithBinning(27);
            _scorer = scorer;

            _seedEnvelopes = new SortedList<double, ObservedIsotopeEnvelope>();
            _seedSet = new List<ObservedIsotopeEnvelope>[NRows];
            for(var i = 0; i < NRows; i++) _seedSet[i] = new List<ObservedIsotopeEnvelope>();
        }

        private LcMsFeatureLikelihood _scorer;        
        public IList<LcMsPeakCluster> FindFeatures(int binNumber)
        {
            var features = GetLcMs1PeakClusters(binNumber);
            return features;
        }

        public Ms1Spectrum GetSpectrum(int ms1ScanNum)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var idx = Array.BinarySearch(ms1ScanNums, ms1ScanNum);
            return idx < 0 ? null : Ms1Spectra[idx];
        }
        
        public LcMsPeakCluster GetLcMsPeakCluster(double targetMass, int targetCharge, int targetMinScanNum, int targetMaxScanNum, bool featureMatrixCreated = false)
        {
            if (!featureMatrixCreated) BuildFeatureMatrix(targetMass); // should be called first

            // 1) set initial elution period for initial sampling
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            if (Run.GetMsLevel(targetMinScanNum) > 1) targetMinScanNum = Run.GetPrevScanNum(targetMinScanNum, 1);
            if (Run.GetMsLevel(targetMaxScanNum) > 1) targetMaxScanNum = Run.GetNextScanNum(targetMaxScanNum, 1);
            
            var row = Math.Max(Math.Min(targetCharge - MinScanCharge, MaxScanCharge - MinScanCharge), 0);
            var minCol = ms1ScanNumToIndex[targetMinScanNum];
            var maxCol = ms1ScanNumToIndex[targetMaxScanNum];
            
            // 2) determine charge state range and elution period
            var range = DetermineFeatureRange(targetMass, row, minCol, maxCol);

            if (range == null) return null;

            var minRow = range.Item1;
            var maxRow = range.Item2;
            minCol = range.Item3;
            maxCol = range.Item4;

            // 3) collect envelopes
          
            var feature = CollectLcMsPeaks(minRow, maxRow, minCol, maxCol);

            if (feature == null) return null;

            feature.UpdateScore(Ms1Spectra);
            if (_scorer != null) feature.Score = _scorer.GetScore(feature);
            else feature.Score = 0d;

            // 4) determine abundance

            var simCutoff = GetDistanceCorrelationThreshold(minRow, maxRow, minCol, maxCol);
            var bcCutoff = Math.Max(0.04, simCutoff.Item1);
            var corrCutoff = Math.Min(0.7, simCutoff.Item2);

            Array.Clear(_xic, 0, _xic.Length);
            foreach (var envelope in feature.EnumerateEnvelopes())
            {
                // skip charge states having no good summed envelope
                //if (_distProfileAcrossCharge[envelope.Charge - MinScanCharge] > GetBcDistThreshold() && _corrProfileAcrossCharge[envelope.Charge - MinScanCharge] < GetCorrThreshold()) continue;
                var envCol = ms1ScanNumToIndex[envelope.ScanNum];
                if (_featureMatrix[envelope.Charge - MinScanCharge][envCol].DivergenceDist > bcCutoff &&
                    _featureMatrix[envelope.Charge - MinScanCharge][envCol].CorrelationCoeff < corrCutoff) continue;
                
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

            // 5) determine abundance only for best charges
            Array.Clear(_xic, 0, _xic.Length);
            foreach (var c in feature.BestCharge)
            {
                var i = c - feature.MinCharge;
                if (i < 0) continue;
                
                foreach (var envelope in feature.Envelopes[i])
                {
                    if (envelope == null) continue;
                    var envCol = ms1ScanNumToIndex[envelope.ScanNum];

                    if (_featureMatrix[envelope.Charge - MinScanCharge][envCol].DivergenceDist > bcCutoff &&
                        _featureMatrix[envelope.Charge - MinScanCharge][envCol].CorrelationCoeff < corrCutoff) continue;

                    _xic[envCol] += envelope.Abundance;
                }
            }
            
            smoothedXic = Smoother.Smooth(_xic);
            abundance = 0d;
            for (var k = 0; k < smoothedXic.Length - 1; k++)
            {
                var centerIntensity = 0.5 * (Math.Max(0, smoothedXic[k]) + Math.Max(0, smoothedXic[k + 1]));
                if (!(centerIntensity > 0)) continue;
                var timeInterval = Run.GetElutionTime(ms1ScanNums[k + 1]) - Run.GetElutionTime(ms1ScanNums[k]);
                abundance += centerIntensity * timeInterval;
            }
            feature.AbundanceForBestCharges = abundance;

            return feature;
        }

        public LcMsRun Run;
        public readonly List<Ms1Spectrum> Ms1Spectra;
        public const int MinScanCharge = 2;
        public const int MaxScanCharge = 60;

        public readonly int NColumns;
        public readonly int NRows;
        public readonly MzComparerWithBinning Comparer;

        public const double RelativeIsotopePeakIntensityThreshold = 0.1d;

        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);

        private readonly int _maxThreadCount;
        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;

        private readonly double[] _distProfileAcrossCharge;
        private readonly double[] _corrProfileAcrossCharge;
        private readonly int[,] _summedEnvelopeColRange;
        private readonly double[] _xic;

        private LcMsPeakMatrixCell[][] _featureMatrix;
        private int[] _rows;
        private int[] _cols;
        private TheoreticalIsotopeEnvelope _theoreticalEnvelope;

        private readonly List<ObservedIsotopeEnvelope>[] _seedSet;
        private readonly SortedList<double, ObservedIsotopeEnvelope> _seedEnvelopes;
        
        private void SetTargetMass(double targetMass)
        {
            var chargeLowerBound = Math.Ceiling(targetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Floor(targetMass / Run.MinMs1Mz);
            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);
            _targetMass = targetMass;
            _rows = Enumerable.Range((int)rowLb, (int)(rowUb - rowLb + 1)).ToArray();
            _theoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, RelativeIsotopePeakIntensityThreshold);
        }
        
        private double _targetMass;

        private void InitFeatureMatrix()
        {
            if (_featureMatrix != null) return;

            _featureMatrix = new LcMsPeakMatrixCell[NRows][];
            for (var i = 0; i < NRows; i++)
            {
                _featureMatrix[i] = new LcMsPeakMatrixCell[NColumns];
                for (var j = 0; j < NColumns; j++)
                    _featureMatrix[i][j] = new LcMsPeakMatrixCell(MaxEnvelopeLength);
            }
        }

        //public const double SNRthreshold = 1.4826;
        private void BuildFeatureMatrix(double targetMass)
        {
            InitFeatureMatrix();

            SetTargetMass(targetMass);
            
            var observedRows        = new bool[NRows];
            var observedCols        = new bool[NColumns];
            
            var mostAbuInternalIdx = _theoreticalEnvelope.IndexOrderByRanking[0];
            var totalElutionLength = Run.GetElutionTime(Run.MaxLcScan);
            var elutionSamplingHalfLen = Math.Max(Math.Min(totalElutionLength * 0.003, 5.0), 0.5);
            var neighborHalfColumns = (int) Math.Max((elutionSamplingHalfLen/totalElutionLength)*NColumns, 5);

            var targetMassBinNum = Comparer.GetBinNumber(targetMass);
            var tolerance = new Tolerance(8);

            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;
            
            var nPeaksCutoff = NumberOfPeaksCutoff;
            var bcSeedCutoff = GetSeedBcDistThreshold();
            var corrSeedCutoff = GetSeedCorrThreshold();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            
            var options = new ParallelOptions();
            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            Parallel.ForEach(_rows, options, row =>
            {
                for (var col = 0; col < NColumns; col++) _featureMatrix[row][col].Init();

                var charge = row + MinScanCharge;

                _seedSet[row].Clear();
               
                for (var k = 0; k < _theoreticalEnvelope.Size; k++)
                {
                    var i = _theoreticalEnvelope.IndexOrderByRanking[k];
                    var isotopeIndex = _theoreticalEnvelope.Isotopes[i].Index;
                    var isotopeMzLb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzStart(targetMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(targetMassBinNum - 1), charge, isotopeIndex);
                    var isotopeMzUb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzEnd(targetMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(targetMassBinNum + 1), charge, isotopeIndex);

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
                            if (_featureMatrix[row][col].EnvelopePeaks[i] == null || ms1Peak.Intensity > _featureMatrix[row][col].EnvelopePeaks[i].Intensity)
                            {
                                _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                _featureMatrix[row][col].AccurateMass = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                            }
                        }
                        else
                        {
                            //if (!(_featureMatrix[row][col].AccurateMass > 0)) _featureMatrix[row][col].AccurateMass = targetMass;
                            if (!(_featureMatrix[row][col].AccurateMass > 0)) continue;
                            var expectedPeakMz = Ion.GetIsotopeMz(_featureMatrix[row][col].AccurateMass, charge, isotopeIndex);
                            if (Math.Abs(expectedPeakMz - ms1Peak.Mz) > tolerance.GetToleranceAsTh(ms1Peak.Mz)) continue;
                            
                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (_featureMatrix[row][col].EnvelopePeaks[i] != null)
                            {
                                if (_featureMatrix[row][col].CountActivePeaks == 1)
                                {
                                    if (ms1Peak.Intensity > _featureMatrix[row][col].EnvelopePeaks[i].Intensity) _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                }
                                else
                                {
                                    var tmpPeak = _featureMatrix[row][col].EnvelopePeaks[i];
                                    var bc1 = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                                    _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                    var bc2 = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                                    if (bc1 < bc2) _featureMatrix[row][col].EnvelopePeaks[i] = tmpPeak;                                    
                                }
                            }
                            else
                            {
                                _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                            }
                        }
                    }

                    if (k == 0)
                    {
                        // for cells missing most abundant peaks
                        for (var col = 0; col < NColumns; col++)
                        {
                            if (_featureMatrix[row][col].Exist) continue;
                            var highestIntensity = 0d;
                            var inferredAccurateMass = 0d;
                            // find the most intense abundant peak from neighboring cells
                            for (var j = Math.Max(col - neighborHalfColumns, 0); j <= Math.Min(col + neighborHalfColumns, NColumns - 1); j++)
                            {
                                var mostAbuPeak = _featureMatrix[row][j].EnvelopePeaks[mostAbuInternalIdx];
                                if (mostAbuPeak != null && mostAbuPeak.Intensity > highestIntensity)
                                {
                                    highestIntensity = mostAbuPeak.Intensity;
                                    inferredAccurateMass = _featureMatrix[row][j].AccurateMass;
                                }
                            }
                            _featureMatrix[row][col].AccurateMass = inferredAccurateMass;
                        }
                    }

                }

                for (var col = 0; col < NColumns; col++)
                {
                    if (!(_featureMatrix[row][col].Exist)) continue;

                    //var mostAbuPeak = _featureMatrix[row][col].EnvelopePeaks[mostAbuIsotopeInternalIndex];
                    //var refIntensity = Ms1Spectra[mostAbuPeak.Ms1SpecIndex].GetLocalMedianIntensity(mostAbuPeak, targetMass);
                    //var signalToNoiseRatio = mostAbuPeak.Intensity / refIntensity;
                    
                    if (_featureMatrix[row][col].CountActivePeaks >= nPeaksCutoff)
                    {
                        var corr = _theoreticalEnvelope.GetPearsonCorrelation(_featureMatrix[row][col].EnvelopePeaks);
                        var bcDist = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                        _featureMatrix[row][col].CorrelationCoeff = corr;
                        _featureMatrix[row][col].DivergenceDist = bcDist;

                        if (!observedRows[row]) observedRows[row] = true;
                        if (!observedCols[col]) observedCols[col] = true;
                        
                        // collect seed envelopes
                        var mostAbuPeak = _featureMatrix[row][col].EnvelopePeaks[mostAbuInternalIdx];
                        if (mostAbuPeak != null && (bcDist < bcSeedCutoff || corr < corrSeedCutoff))
                        {
                            var signalToNoiseRatio = mostAbuPeak.Intensity / Ms1Spectra[col].MedianIntensity;
                            if (signalToNoiseRatio > 3)
                            {
                                var seed = new ObservedIsotopeEnvelope(_featureMatrix[row][col].AccurateMass, row + MinScanCharge, ms1ScanNums[col], _featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope);
                                _seedSet[row].Add(seed);
                            }
                        }
                    }
                    else
                    {
                        _featureMatrix[row][col].AccurateMass = 0d;
                    }
                }
            }// end or row for-loop
            );

            var temp = new List<int>();
            for (var i = 0; i < observedRows.Length; i++) if (observedRows[i]) temp.Add(i);
            _rows = temp.ToArray();

            temp.Clear();
            for (var i = 0; i < observedCols.Length; i++) if (observedCols[i]) temp.Add(i);
            _cols = temp.ToArray();


            // sort seed envelopes
            _seedEnvelopes.Clear();
            foreach (var row in _rows)
            {
                foreach (var seed in _seedSet[row])
                {
                    var col = ms1ScanNumToIndex[seed.ScanNum];
                    var bcDist = _featureMatrix[row][col].DivergenceDist;
                    _seedEnvelopes.Add(bcDist, seed);
                }
            }
        }

        private int NumberOfPeaksCutoff
        {
            get
            {
                if (_targetMass < 2000) return 2;
                if (_targetMass < 8000) return 3;
                if (_targetMass < 25000) return 4;
                if (_targetMass < 50000) return 5;
                return 6;
            }
        }
        /*
        private IList<ObservedIsotopeEnvelope> GetSeedCells()
        {
            var seedList = new List<KeyValuePair<double, ObservedIsotopeEnvelope>>();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var bcCutoff = GetSeedBcDistThreshold();
            var corrCutoff = GetSeedCorrThreshold();

            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            
            foreach (var i in _rows)
            {
                foreach (var j in _cols)
                {
                    if (!_featureMatrix[i][j].Exist) continue;
                    
                    var mostAbuPeak = _featureMatrix[i][j].EnvelopePeaks[mostAbuInternalIndex];
                    if (mostAbuPeak == null) continue;

                    if (_featureMatrix[i][j].CountActivePeaks < NumberOfPeaksCutoff) continue;

                    var bcDist = _featureMatrix[i][j].DivergenceDist;
                    var corr = _featureMatrix[i][j].CorrelationCoeff;
                    
                    if (bcDist > bcCutoff && corr < corrCutoff) continue;

                    var signalToNoiseRatio = mostAbuPeak.Intensity / Ms1Spectra[j].MedianIntensity;
                    if (signalToNoiseRatio < 3) continue;
                    
                    var seed = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, _theoreticalEnvelope);
                    seedList.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(bcDist, seed));
                    //seedList.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(corr, seed));
                }
            }
            //return seedList.OrderByDescending(x => x.Key).Select(x => x.Value).ToList();
            return seedList.OrderBy(x => x.Key).Select(x => x.Value).ToList();
        }*/

        private IList<LcMsPeakCluster> GetLcMs1PeakClusters(int binNumber)
        {
            const int chargeNeighborGap = 4;
            var targetMass = Comparer.GetMzAverage(binNumber);
            BuildFeatureMatrix(targetMass); // should be called first

            var clusters = new List<LcMsPeakCluster>();

            // todo : bottom up dataset??
            if (_rows.Length < 2) return clusters;

            var tempEnvelope = new double[_theoreticalEnvelope.Size];
            var tempEnvelope2 = new double[_theoreticalEnvelope.Size];
            
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            var tolerance = new Tolerance(5);
            //var seedList = GetSeedCells();

            foreach (var seedItem in _seedEnvelopes)
            {
                var seed = seedItem.Value;
                var row = seed.Charge - MinScanCharge;
                var col = ms1ScanNumToIndex[seed.ScanNum];
                
                if (_featureMatrix[row][col].CheckedOutFlag) continue;
                
                var mostAbuMz = _theoreticalEnvelope.GetIsotopeMz(row + MinScanCharge, mostAbuInternalIndex);
                var seedLocalWin = Ms1Spectra[col].GetLocalMzWindow(mostAbuMz);
                var poissonPvalue = seedLocalWin.GetPoissonTestPvalue(_featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope.Size);
                var rankSumPvalue = seedLocalWin.GetRankSumTestPvalue(_featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope.Size);

                //var goodEnvelope = (rankSumPvalue < 0.01 && poissonPvalue < 0.01) || (rankSumPvalue < 1e-3) || (poissonPvalue < 1e-3);
                var goodEnvelope = (rankSumPvalue < 0.01 || poissonPvalue < 0.01);
                if (!goodEnvelope) continue;

                var chargeCheck = CorrectChargeState(seed, Ms1Spectra[col]);
                if (!chargeCheck) continue;
                
                var seedMass = _featureMatrix[row][col].AccurateMass;
                var massTol = tolerance.GetToleranceAsTh(seedMass);
                var newCluster = new LcMsPeakCluster(Run, seed);

                Array.Clear(tempEnvelope, 0, tempEnvelope.Length);
                seed.Peaks.SumEnvelopeTo(tempEnvelope);

                var neighbors = new Queue<ObservedIsotopeEnvelope>();
                neighbors.Enqueue(seed); // pick a seed
                _featureMatrix[row][col].CheckedOutFlag = true;

                var summedBcDist = _featureMatrix[row][col].DivergenceDist;
                var summedCorr = _featureMatrix[row][col].CorrelationCoeff;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    var charge = cell.Charge;
                    
                    var minRw = (int) Math.Max(charge - MinScanCharge - chargeNeighborGap, _rows.First());
                    var maxRw = (int) Math.Min(charge - MinScanCharge + chargeNeighborGap, _rows.Last());
                    var currCol = ms1ScanNumToIndex[cell.ScanNum];

                    for (var k = 0; k < 5; k++)
                    {
                        var j = currCol;
                        if (k < 3) j += k;
                        else j -= (k - 2);

                        if (j < _cols.First() || j > _cols.Last()) continue;

                        for (var i = minRw; i <= maxRw; i++)
                        {
                            if (_featureMatrix[i][j].CheckedOutFlag) continue;
                            if (!(_featureMatrix[i][j].AccurateMass > 0)) continue;
                            if (Math.Abs(seedMass - _featureMatrix[i][j].AccurateMass) > massTol) continue;
                            
                            Array.Copy(tempEnvelope, tempEnvelope2, tempEnvelope2.Length);

                            _featureMatrix[i][j].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                            var newDivergence = _theoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                            var newCorrelation = _theoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);
                            
                            if (_featureMatrix[i][j].DivergenceDist < 0.02 &&
                                _featureMatrix[i][j].CorrelationCoeff > 0.8 || newDivergence < summedBcDist ||
                                newCorrelation > summedCorr)
                            {
                                var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass,
                                    i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks,
                                    _theoreticalEnvelope);

                                neighbors.Enqueue(envelope);
                                newCluster.Expand(envelope);
                                _featureMatrix[i][j].CheckedOutFlag = true;

                                summedBcDist = newDivergence;
                                summedCorr = newCorrelation;
                            }
                            else
                            {
                                Array.Copy(tempEnvelope2, tempEnvelope, tempEnvelope.Length);
                            }
                        }
                    }
                }

                LcMsPeakCluster refinedCluster = null;
                if (summedCorr > 0.5 || summedBcDist < 0.15)
                {
                    // re-update check-out map
                    SetCheckOutFlag(newCluster.MinCharge - MinScanCharge, newCluster.MaxCharge - MinScanCharge, ms1ScanNumToIndex[newCluster.MinScanNum], ms1ScanNumToIndex[newCluster.MaxScanNum], false);
                    refinedCluster = GetLcMsPeakCluster(newCluster.RepresentativeMass, newCluster.RepresentativeCharge, newCluster.MinScanNum, newCluster.MaxScanNum, true);
                }

                if (refinedCluster != null && refinedCluster.GoodEnougth && (_scorer != null && refinedCluster.Score >= _scorer.ScoreThreshold))
                {
                    SetCheckOutFlag(_rows.First(), _rows.Last(), ms1ScanNumToIndex[refinedCluster.MinScanNum], ms1ScanNumToIndex[refinedCluster.MaxScanNum], true);
                    clusters.Add(refinedCluster);
                }
                else
                {
                    SetCheckOutFlag(newCluster.MinCharge - MinScanCharge, newCluster.MaxCharge - MinScanCharge, ms1ScanNumToIndex[newCluster.MinScanNum], ms1ScanNumToIndex[newCluster.MaxScanNum], true);
                }
            }
            return clusters;            
        }

        private void SetCheckOutFlag(int minRow, int maxRow, int minCol, int maxCol, bool flag)
        {
            for(var i = minRow; i <= maxRow; i++)
                for (var j = minCol; j <= maxCol; j++) _featureMatrix[i][j].CheckedOutFlag = flag;
        }

        private Tuple<double, double> GetDistanceCorrelationThreshold(int minRow, int maxRow, int minCol, int maxCol)
        {
            var distList = new List<double>();
            var corrList = new List<double>();

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist) continue;
                    distList.Add(_featureMatrix[i][j].DivergenceDist);
                    corrList.Add(_featureMatrix[i][j].CorrelationCoeff);
                }
            }

            //var meanStd = distList.MeanStandardDeviation();
            //var bcCutoff = Math.Min(Math.Max(meanStd.Item1 - meanStd.Item2 * sigma, 0.05), 0.3);
            var distCutoff = distList.Median();
            var corrCutoff = corrList.Median();

            return new Tuple<double, double>(distCutoff, corrCutoff);
        }
        
        private LcMsPeakCluster CollectLcMsPeaks(int minRow, int maxRow, int minCol, int maxCol)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var envelopes = new List<ObservedIsotopeEnvelope>();
            var bestBcDist = 100d;
            ObservedIsotopeEnvelope bestEnvelope = null;
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            var bcCutoff = GetSeedBcDistThreshold();
            var corrCutoff = GetSeedCorrThreshold();
            //var simCutoff = GetDistanceCorrelationThreshold(minRow, maxRow, minCol, maxCol);
            //var bcCutoff = Math.Max(0.04, simCutoff.Item1);
            //var corrCutoff = Math.Min(0.7, simCutoff.Item2);

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist) continue;

                    if (_featureMatrix[i][j].DivergenceDist > bcCutoff && _featureMatrix[i][j].CorrelationCoeff < corrCutoff) continue; // exclude outliers

                    var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, _theoreticalEnvelope);
                    envelopes.Add(envelope);

                    if (_featureMatrix[i][j].EnvelopePeaks[mostAbuInternalIndex] != null && _featureMatrix[i][j].DivergenceDist < bestBcDist)
                    {
                        bestBcDist = _featureMatrix[i][j].DivergenceDist;
                        bestEnvelope = envelope;
                    }
                }
            }

            if (bestEnvelope == null) return null;

            var cluster = new LcMsPeakCluster(Run, bestEnvelope);
            cluster.AddEnvelopes(minRow + MinScanCharge, maxRow + MinScanCharge, ms1ScanNums[minCol], ms1ScanNums[maxCol], envelopes);
            //foreach(var e in envelopes) if (bestEnvelope != e) cluster.AddObservedEnvelope(e);
            
            return cluster;
        }


        private double GetBcDistThreshold()
        {
            if (_targetMass > 35000) return 0.12;
            if (_targetMass > 25000) return 0.1;
            if (_targetMass > 15000) return 0.08;
            return 0.06;

        }

        private double GetCorrThreshold()
        {
            if (_targetMass > 35000) return 0.3;
            if (_targetMass > 25000) return 0.4;
            if (_targetMass > 15000) return 0.5;
            return 0.6;
        }


        private double GetSeedBcDistThreshold()
        {
            if (_targetMass > 45000) return 0.3;
            if (_targetMass > 35000) return 0.25;
            if (_targetMass > 25000) return 0.2;
            if (_targetMass > 15000) return 0.15;
            return 0.1;
        }

        private double GetSeedCorrThreshold()
        {
            if (_targetMass > 45000) return 0.2;
            if (_targetMass > 35000) return 0.3;
            if (_targetMass > 25000) return 0.4;
            if (_targetMass > 15000) return 0.5;
            return 0.6;
        }

        private Tuple<int, int, int, int> DetermineFeatureRange(double targetMass, int seedRow, int seedMinCol, int seedMaxCol)
        {
            Array.Clear(_distProfileAcrossCharge, 0, _distProfileAcrossCharge.Length);
            Array.Clear(_corrProfileAcrossCharge, 0, _corrProfileAcrossCharge.Length);

            // expand seed min/max col when the period is too short            
            var elutionSamplingPeriod = Math.Max(Math.Min(Run.GetElutionTime(Run.MaxLcScan) * 0.003, 5.0), 0.5);
            var tempWindow = GetElutionWindow((int)(0.5*(seedMinCol+seedMaxCol)), elutionSamplingPeriod);
            for (var t = seedMinCol - 1; t >= tempWindow.Item1; t--)
            {
                if (_featureMatrix[seedRow][t].CheckedOutFlag) break;
                seedMinCol = t;
            }
            for (var t = seedMaxCol + 1; t <= tempWindow.Item2; t++)
            {
                if (_featureMatrix[seedRow][t].CheckedOutFlag) break;
                seedMaxCol = t;
            }

            var bcCutoff = GetBcDistThreshold();
            var corrCutoff = GetCorrThreshold();
            var bcDistList = new List<double>();
            var corrList = new List<double>();


            var options = new ParallelOptions();
            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            //foreach (var i in _rows)
            Parallel.ForEach(_rows, options, i =>
            {
                var summedEnvelope = GetSummedEnvelopeAtCharge(targetMass, i, seedMinCol, seedMaxCol);
                //if (_corrProfileAcrossCharge[i] < corrCutoff && _distProfileAcrossCharge[i] > bcCutoff) continue;
                //bcDistList.Add(_distProfileAcrossCharge[i]);
                //corrList.Add(_corrProfileAcrossCharge[i]);
            });
            
            foreach (var i in _rows)
            {
                //var summedEnvelope = GetSummedEnvelopeAtCharge(targetMass, i, seedMinCol, seedMaxCol);
                if (_corrProfileAcrossCharge[i] < corrCutoff && _distProfileAcrossCharge[i] > bcCutoff) continue;
                bcDistList.Add(_distProfileAcrossCharge[i]);
                corrList.Add(_corrProfileAcrossCharge[i]);                
            }

            var minCol = seedMinCol;
            var maxCol = seedMaxCol;
            if (bcDistList.Count > 0)
            {
                bcDistList.Sort();
                corrList.Sort();
                var bcDistMedian = bcDistList[(int)(bcDistList.Count * 0.5)];
                var corrMedian = corrList[(int)(corrList.Count * 0.5)];
                var minColList = new List<int>();
                var maxColList = new List<int>();
                foreach (var i in _rows)
                {
                    if (_distProfileAcrossCharge[i] > bcDistMedian && _corrProfileAcrossCharge[i] < corrMedian) continue;
                    minColList.Add(_summedEnvelopeColRange[i, 0]);
                    maxColList.Add(_summedEnvelopeColRange[i, 1]);
                }
                minColList.Sort();
                maxColList.Sort();
                minCol = minColList[(int) (minColList.Count * 0.5)];
                maxCol = maxColList[(int) (maxColList.Count * 0.5)];
            }

            // for each carge state, find and evaluate optimial summed envelope
            // find column range
            var bestCharge = new int[] { 0, 0 };
            var bestChargeDist = new double[] { 0, 0 };
            var minRow = _rows.Last();
            var maxRow = _rows.First();
            foreach (var i in _rows)
            {
                if (_corrProfileAcrossCharge[i] < corrCutoff && _distProfileAcrossCharge[i] > bcCutoff) continue;

                if (_summedEnvelopeColRange[i, 1] < minCol) continue;
                if (_summedEnvelopeColRange[i, 0] > maxCol) continue;

                if (i < minRow) minRow = i;
                if (i > maxRow) maxRow = i;

                var charge = i + MinScanCharge;
                var chargeIdx = (charge % 2 == 0) ? 0 : 1;
                if (bestCharge[chargeIdx] == 0 || _distProfileAcrossCharge[i] < bestChargeDist[chargeIdx])
                {
                    bestChargeDist[chargeIdx] = _distProfileAcrossCharge[i];
                    bestCharge[chargeIdx] = charge;
                }
                //if (_summedEnvelopeColRange[i,0] < minCol) minCol = _summedEnvelopeColRange[i,0];
                //if (_summedEnvelopeColRange[i,1] > maxCol) maxCol = _summedEnvelopeColRange[i,1];
            }
            
            // only one charge? there should be another....force to cover neighboring charge states
            if (bestCharge[0] == 0 && bestCharge[1] == 0) return null;

            if (bestCharge[0] == 0 || bestCharge[1] == 0)
            {
                minRow = Math.Max(minRow - 1, _rows.First());
                maxRow = Math.Min(maxRow + 1, _rows.Last());
            }

            if (minRow == maxRow) return null;

            for (var i = minRow; i <= maxRow; i++)
            {
                var charge = i + MinScanCharge;
                var chargeIdx = (charge % 2 == 0) ? 0 : 1;
                if (bestCharge[chargeIdx] == 0)
                {
                    bestCharge[chargeIdx] = charge;
                    break;
                }
            }

            // refine scan boundary by considering bell shaped Xic for the best charge state
            Array.Clear(_xic, 0, _xic.Length);
            var nIntensity = 0;
            foreach (var c in bestCharge)
            {
                var i = c - MinScanCharge;
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist ||
                        (_featureMatrix[i][j].DivergenceDist > GetSeedBcDistThreshold() && _featureMatrix[i][j].CorrelationCoeff < GetSeedCorrThreshold())) continue;
                    _xic[j] += _featureMatrix[i][j].Intensity;
                    nIntensity++;
                }
            }
            /*
            var bestBcDistRow = (bestChargeAbundance[0] > bestChargeAbundance[1]) ? bestCharge[0] - MinScanCharge : bestCharge[1] - MinScanCharge;
            for (var j = minCol; j <= maxCol; j++)
            {
                if (!_featureMatrix[bestBcDistRow][j].Exist || 
                    (_featureMatrix[bestBcDistRow][j].DivergenceDist > GetSeedBcDistThreshold() && _featureMatrix[bestBcDistRow][j].CorrelationCoeff < GetSeedCorrThreshold())) continue;

                _xic[j] = _featureMatrix[bestBcDistRow][j].Intensity;
                nIntensity++;
            } */

            if (nIntensity < 2)
            {
                //return new Tuple<int, int, int, int>(minRow, maxRow, oriSeedMinCol, oriSeedMaxCol);
                return null;
            }

            var smoothedXic = Smoother.Smooth(_xic);
            var apexCol = -1;
            var apexIntensity = 0d;
            for (var j = minCol; j <= maxCol; j++)
            {
                if (smoothedXic[j] > apexIntensity)
                {
                    apexIntensity = smoothedXic[j];
                    apexCol = j;
                }
            }

            if (apexCol < 0) return null;
            
            ////////////////// 1) First Half
            var ms1ScanNums = Run.GetMs1ScanVector();
            var oneSigIntensity = apexIntensity * Math.Exp(-1);
            var threeSigIntensity = apexIntensity * Math.Exp(-4.5);
            
            // estimate sigma for Gaussian shaped elution profile for the first half
            var elutionStartColByOneSigma = apexCol;
            for (var j = apexCol - 1; j >= 0; j--)
            {
                elutionStartColByOneSigma = j;
                if (smoothedXic[j] < oneSigIntensity) break;
            }

            var oneSigPeriod = Run.GetElutionTime(ms1ScanNums[apexCol]) -
                               Run.GetElutionTime(ms1ScanNums[elutionStartColByOneSigma]);

            // find 2-sigma position at elution time prior to Apex (going backward)
            var elutionStartColByTwoSigma = elutionStartColByOneSigma;
            for (var j = elutionStartColByOneSigma - 1; j >= 0; j--)
            {
                elutionStartColByTwoSigma = j;
                if (Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]) > 2*oneSigPeriod) break;
            }

            // extends for long tail
            var elutionStartCol = elutionStartColByTwoSigma;
            //for (var j = elutionStartColByTwoSigma - 1; j >= 0; j--)
            for (var j = elutionStartColByTwoSigma - 1; j >= minCol; j--)
            {
                if (smoothedXic[j] < threeSigIntensity) break;
                if (Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]) > 3 * oneSigPeriod) break;

                var needBreak = true;
                for (var k = j; k >= j - 3; k--)
                {
                    if (k < 0) break;
                    if (smoothedXic[k] < smoothedXic[k + 1]) needBreak = false;
                }
                if (needBreak) break;
                    
                elutionStartCol = j;
            }

            ////////////////// Second Half
            var elutionEndColByOneSigma = apexCol;
            for (var j = apexCol + 1; j < NColumns; j++)
            {
                elutionEndColByOneSigma = j;
                if (smoothedXic[j] < oneSigIntensity) break;
            }
            // find two sigma point at elution time axis
            oneSigPeriod = Run.GetElutionTime(ms1ScanNums[elutionEndColByOneSigma]) - Run.GetElutionTime(ms1ScanNums[apexCol]);
            var elutionEndColByTwoSigma = elutionEndColByOneSigma;
            for (var j = elutionEndColByOneSigma + 1; j < NColumns; j++)
            {
                elutionEndColByTwoSigma = j;
                if (Run.GetElutionTime(ms1ScanNums[j]) -  Run.GetElutionTime(ms1ScanNums[apexCol]) > 2 * oneSigPeriod) break;
            }

            // possible extension?
            var elutionEndCol = elutionEndColByTwoSigma;
            //for (var j = elutionEndColByTwoSigma + 1; j < NColumns; j++)
            for (var j = elutionEndColByTwoSigma + 1; j <= maxCol; j++)
            {
                if (smoothedXic[j] < threeSigIntensity) break;
                if (Run.GetElutionTime(ms1ScanNums[j]) - Run.GetElutionTime(ms1ScanNums[apexCol]) > 4 * oneSigPeriod) break;

                var needBreak = true;
                for (var k = j; k <= j + 3; k++)
                {
                    if (k >= NColumns) break;
                    if (smoothedXic[k] < smoothedXic[k - 1]) needBreak = false;
                }
                if (needBreak) break;

                elutionEndCol = j;
            }
            
            return new Tuple<int, int, int, int>(minRow, maxRow, elutionStartCol, elutionEndCol);
        }

        private double[] GetSummedEnvelopeAtCharge(double targetMass, int row, int minCol, int maxCol)
        {
            const int maxScanSkips = 2;
            const double goodEnoughBcDistance = 0.04;
            const double goodEnoughCorrCoeff = 0.7;
            var summedEnvelope = new double[_theoreticalEnvelope.Size];
            var seedCol = -1;
            var seedDist = 10.0d;
            
            var newMinCol = minCol;
            var newMaxCol = maxCol;

            var tolerance = new Tolerance(5);
            var massTol = tolerance.GetToleranceAsTh(targetMass);

            for (var j = minCol; j <= maxCol; j++)
            {
                if (!_featureMatrix[row][j].Exist) continue;
                if (_featureMatrix[row][j].CheckedOutFlag) continue;
                if (Math.Abs(targetMass - _featureMatrix[row][j].AccurateMass) > massTol) continue;
                
                var bcDist = _featureMatrix[row][j].DivergenceDist;
                if (bcDist > seedDist) continue;

                //var signalToNoiseRatio = _featureMatrix[row][j].HighestIntensity / Ms1Spectra[j].MedianIntensity;
                //if (signalToNoiseRatio < 3) continue;
                seedCol = j;
                seedDist = bcDist;
                
            }

            var summedBcDist = 1.0d;
            var summedCorr = 0d;

            _corrProfileAcrossCharge[row] = 0d;
            _distProfileAcrossCharge[row] = 1.0d;

            if (seedCol < 0) return summedEnvelope;

            // going forward
            var tempEnvelope = new double[_theoreticalEnvelope.Size];
            
            var n = 0;
            for(var col = seedCol; col < NColumns; col++)
            {
                if (_featureMatrix[row][col].CheckedOutFlag) break;

                if (n >= maxScanSkips)
                {
                    newMaxCol = col;
                    break;
                }

                if (!_featureMatrix[row][col].Exist)
                {
                    n++;
                    continue;
                }

                _featureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                var tempBcDist = _theoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                var tempCorr = _theoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);

                if (tempBcDist < summedBcDist || tempCorr > summedCorr || _featureMatrix[row][col].DivergenceDist < goodEnoughBcDistance || _featureMatrix[row][col].CorrelationCoeff > goodEnoughCorrCoeff)
                {
                    summedBcDist = tempBcDist;
                    summedCorr = tempCorr;
                    Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                    n = 0;
                }
                else
                {
                    Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                    n++;
                }
            }

            // going backward
            summedBcDist = 10.0d;
            summedCorr = 0d;
            Array.Clear(tempEnvelope, 0, tempEnvelope.Length);
            Array.Clear(summedEnvelope, 0, summedEnvelope.Length);
            n = 0;
            for (var col = seedCol; col >= 0; col--)
            {
                if (_featureMatrix[row][col].CheckedOutFlag) break;

                if (n >= maxScanSkips)
                {
                    newMinCol = col;
                    break;
                }

                if (!_featureMatrix[row][col].Exist)
                {
                    n++;
                    continue;
                }

                _featureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                var tempBcDist = _theoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                var tempCorr = _theoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);

                if (tempBcDist < summedBcDist || tempCorr > summedCorr || _featureMatrix[row][col].DivergenceDist < goodEnoughBcDistance || _featureMatrix[row][col].CorrelationCoeff > goodEnoughCorrCoeff)
                {
                    summedBcDist = tempBcDist;
                    summedCorr = tempCorr;
                    Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                    n = 0;
                }
                else
                {
                    Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                    n++;
                }
            }

            // construct summed envelope, given newMinCol and newMaxCol
            summedBcDist = 10.0d;
            summedCorr = 0d;
            Array.Clear(tempEnvelope, 0, tempEnvelope.Length);
            _featureMatrix[row][seedCol].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
            Array.Copy(tempEnvelope, summedEnvelope, summedEnvelope.Length);
            var colShift = 0;
            var hitMinCol = false;
            var hitMaxCol = false;
            while (true)
            {
                if (hitMinCol && hitMaxCol) break;
                if (n > 3) break;
                colShift++;

                for (var colDir = 0; colDir < 2; colDir++)
                {
                    var col = (colDir == 0) ? seedCol + colShift : seedCol - colShift;

                    if (col < newMinCol)
                    {
                        hitMinCol = true;
                        continue;
                    }
                    if (col > newMaxCol)
                    {
                        hitMaxCol = true;
                        continue;
                    }

                    if (!_featureMatrix[row][col].Exist) continue;

                    _featureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                    var tempBcDist = _theoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                    var tempCorr = _theoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);

                    if (tempBcDist < summedBcDist || tempCorr > summedCorr)
                    {
                        summedBcDist = tempBcDist;
                        summedCorr = tempCorr;
                        Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                    }
                    else
                    {
                        Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                        n++;
                    }                    
                }
            }
            
            _corrProfileAcrossCharge[row] = summedCorr;
            _distProfileAcrossCharge[row] = summedBcDist;

            _summedEnvelopeColRange[row, 0] = newMinCol;
            _summedEnvelopeColRange[row, 1] = newMaxCol;
            
            return summedEnvelope;
        }


        private bool CorrectChargeState(ObservedIsotopeEnvelope envelope, Ms1Spectrum spectrum)
        {
            //var checkCharge = envelope.Charge;
            if (envelope.Charge > 20) return true; //high charge (> +20), just pass

            var peaks = spectrum.Peaks;
            var peakStartIndex = envelope.MinMzPeak.IndexInSpectrum;
            var peakEndIndex = envelope.MaxMzPeak.IndexInSpectrum;
            var intensityThreshold = envelope.HighestIntensity * 0.15;
            
            var nPeaks = 0;
            for (var i = peakStartIndex; i <= peakEndIndex; i++)
            {
                if (peaks[i].Intensity > intensityThreshold) nPeaks++;
            }

            //if (nPeaks < 10) return false;
            if (envelope.NumberOfPeaks > nPeaks * 0.7) return true;

            var tolerance = new Tolerance(5);
            var threshold = nPeaks * 0.5;
            var threshold2 = envelope.NumberOfPeaks + (envelope.TheoreticalEnvelope.Size - 1) * 0.7;

            var mzTol = tolerance.GetToleranceAsTh(peaks[peakStartIndex].Mz);

            var minCheckCharge = Math.Max(envelope.Charge * 2 - 1, 4);
            var maxCheckCharge = Math.Min(envelope.Charge * 5 + 1, 60);
            var maxDeltaMz = Constants.C13MinusC12 / minCheckCharge + mzTol;
            var nChargeGaps = new int[maxCheckCharge - minCheckCharge + 1];
            
            for (var i = peakStartIndex; i <= peakEndIndex; i++)
            {
                if (!(peaks[i].Intensity > intensityThreshold)) continue;
                
                for (var j = i + 1; j <= peakEndIndex; j++)
                {
                    if (!(peaks[j].Intensity > intensityThreshold)) continue;
                    
                    var deltaMz = peaks[j].Mz - peaks[i].Mz;

                    if (deltaMz > maxDeltaMz) break;
                    for (var c = Math.Round(1 / (deltaMz + mzTol)); c <= Math.Round(1 / (deltaMz - mzTol)); c++)
                    {
                        if (c < minCheckCharge || c > maxCheckCharge) continue;
                        var k = (int)c - minCheckCharge;
                        nChargeGaps[k]++;

                        //if (nChargeGaps[k] + 1 > threshold && nChargeGaps[k] + 1 > 1.25 * envelope.NumberOfPeaks) return false;
                        if (nChargeGaps[k] + 1 > threshold && nChargeGaps[k] + 1 > threshold2) return false;
                    }
                }
            }

            return true;
        }
        
        private Tuple<int, int> GetElutionWindow(int refCol, double halfWindowSize, int minScanCount = 5)
        {
            var maxCol = (int)Math.Min(minScanCount + refCol, _cols.Last());
            var minCol = (int)Math.Max(refCol - minScanCount, _cols.First());
            maxCol = Math.Max(maxCol, ShiftColoumByElutionTime(refCol, halfWindowSize));
            minCol = Math.Max(minCol, ShiftColoumByElutionTime(refCol, -halfWindowSize));
            return new Tuple<int, int>(minCol, maxCol);
        }

        private int ShiftColoumByElutionTime(int refCol, double elutionShift)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var elutionTime = Run.GetElutionTime(ms1ScanNums[refCol]);

            if (elutionShift > 0)
            {
                for (var j = refCol + 1; j < _cols.Last(); j++)
                {
                    if (Run.GetElutionTime(ms1ScanNums[j]) > elutionTime + elutionShift) return j;
                }
                return _cols.Last();
            }
            else
            {
                for (var j = refCol - 1; j >= _cols.First(); j--)
                {
                    if (Run.GetElutionTime(ms1ScanNums[j]) < elutionTime + elutionShift) return j;
                }
                return _cols.First();
            }
        }
        
        internal class LcMsPeakMatrixCell
        {
            internal double CorrelationCoeff;
            internal double DivergenceDist;
            internal double AccurateMass;

            internal readonly Ms1Peak[] EnvelopePeaks;
            internal bool CheckedOutFlag;

            internal LcMsPeakMatrixCell(int maxPeaks)
            {
                EnvelopePeaks = new Ms1Peak[maxPeaks];
                CorrelationCoeff = 0;
                DivergenceDist = 0;
                AccurateMass = 0;
                CheckedOutFlag = false;
            }

            internal void Init()
            {
                CorrelationCoeff = 0;
                DivergenceDist = 0;
                AccurateMass = 0;
                Array.Clear(EnvelopePeaks, 0, EnvelopePeaks.Length);
                CheckedOutFlag = false;
            }

            internal bool Exist
            {
                get { return (AccurateMass > 0d); }
            }

            internal int CountActivePeaks
            {
                get { return EnvelopePeaks.Count(p => p != null && p.Active); }
            }
            
            internal double Intensity
            {
                get { return EnvelopePeaks.Where(envelopePeak => envelopePeak != null && envelopePeak.Active).Sum(envelopePeak => envelopePeak.Intensity); }
            }

            internal double HighestIntensity
            {
                get { return EnvelopePeaks.Where(peak => peak != null && peak.Active).Max(peak => peak.Intensity); }
            }
        }
    }
}
