using System;
using System.Collections;
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
        public LcMsPeakMatrix(LcMsRun run, ISequenceFilter ms1Filter, LcMsFeatureLikelihood scorer = null, int minScanCharge = 1, int maxScanCharge = 60,
            int maxThreadCount = 0)
            : this(run, scorer, minScanCharge, maxScanCharge, maxThreadCount)
        {
            _ms1Filter = ms1Filter;
        }
        
        public LcMsPeakMatrix(LcMsRun run, LcMsFeatureLikelihood scorer = null, int minScanCharge = 1, int maxScanCharge = 60, int maxThreadCount = 0)
        {
            Run = run;
            _maxThreadCount = maxThreadCount;

            _ms1PeakList = new List<Ms1Peak>();
            Ms1Spectra = new List<Ms1Spectrum>();
            var ms1ScanNums = run.GetMs1ScanVector();

            NColumns = ms1ScanNums.Length;
            MinSearchCharge = minScanCharge;
            MaxSearchCharge = maxScanCharge;
            NRows = Math.Min(MaxSearchCharge - MinSearchCharge + 1, 35);

            for (var i = 0; i < Math.Min(ms1ScanNums.Length, ushort.MaxValue); i++)
            {
                var ms1Spec = run.GetMs1Spectrum(ms1ScanNums[i]);
                Ms1Spectra.Add(ms1Spec);
                _ms1PeakList.AddRange((Ms1Peak[])ms1Spec.Peaks);
            }
            _ms1PeakList.Sort();

            _distProfileAcrossCharge = new double[NRows];
            _corrProfileAcrossCharge = new double[NRows];
            _intensityAcrossCharge = new double[NRows];
            _summedEnvelopeColRange = new int[NRows, 2];
            
            _featureMatrix = null;
            Comparer = new MzComparerWithBinning(27); // 16ppm
            _scorer = scorer;
            _seedEnvelopes = new List<KeyValuePair<double, ObservedIsotopeEnvelope>>();
            _ms1Filter = null;
        }

        private ISequenceFilter _ms1Filter;
        private readonly LcMsFeatureLikelihood _scorer;        
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

            if (_rows.Length < 2 || _cols.Length < 1) return null;

            // 1) set initial elution period for initial sampling
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            if (Run.GetMsLevel(targetMinScanNum) > 1) targetMinScanNum = Run.GetPrevScanNum(targetMinScanNum, 1);
            if (Run.GetMsLevel(targetMaxScanNum) > 1) targetMaxScanNum = Run.GetNextScanNum(targetMaxScanNum, 1);
            
            //var row = Math.Max(Math.Min(targetCharge - MinScanCharge, MaxScanCharge - MinScanCharge), 0);
            var row = targetCharge - _targetMinCharge;
            if (row < _rows.First() || row > _rows.Last()) row = _rows[(int) (_rows.Length*0.5)];

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
            var feature = CollectLcMsPeaks(targetMass, minRow, maxRow, minCol, maxCol);

            if (feature == null) return null;

            // 4) scoring
            feature.UpdateScore(Ms1Spectra);
            if (_scorer != null) feature.Score = _scorer.GetScore(feature);
            else feature.Score = 0d;
           
            // 5) determine abundance
            SetAbundanceByAUC(ref feature);
           
            return feature;
        }

        public LcMsPeakCluster GetLcMsPeakCluster(double targetMass, int minCharge, int maxCharge, int minScanNum, int maxScanNum)
        {
            InitFeatureMatrix();
            SetTargetMass(targetMass);
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            var minCol = ms1ScanNumToIndex[minScanNum];
            var maxCol = ms1ScanNumToIndex[maxScanNum];
            var minRow = minCharge - _targetMinCharge;
            var maxRow = maxCharge - _targetMinCharge;

            // 3) collect envelopes
            var feature = CollectLcMsPeaks(targetMass, minRow, maxRow, minCol, maxCol, true);

            if (feature == null) return null;

            // 4) scoring
            feature.UpdateScore(Ms1Spectra);
            if (_scorer != null) feature.Score = _scorer.GetScore(feature);
            else feature.Score = 0d;

            // 5) determine abundance
            SetAbundanceByAUC(ref feature);

            return feature;
        }
        
        public void SetAbundanceByAUC(ref LcMsPeakCluster feature)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var ms1ScanNums = Run.GetMs1ScanVector();

            var minCol = ms1ScanNumToIndex[feature.MinScanNum];
            var maxCol = ms1ScanNumToIndex[feature.MaxScanNum];
            var minRow = feature.MinCharge - _targetMinCharge;
            var maxRow = feature.MaxCharge - _targetMinCharge;

            var simCutoff = GetDistanceCorrelationThreshold(minRow, maxRow, minCol, maxCol);
            var bcCutoff = Math.Max(0.07, simCutoff.Item1);
            var corrCutoff = Math.Min(0.7, simCutoff.Item2);

            //Array.Clear(_xic, 0, _xic.Length);
            var xic = new double[maxCol - minCol + 1 + 18];
            const int xicStartIndex = 9;
            var xicEndIndex = xic.Length - 9 - 1;

            foreach (var envelope in feature.EnumerateEnvelopes())
            {
                var envCol = ms1ScanNumToIndex[envelope.ScanNum];
                if (_featureMatrix[envelope.Charge - _targetMinCharge][envCol].DivergenceDist > bcCutoff &&
                    _featureMatrix[envelope.Charge - _targetMinCharge][envCol].CorrelationCoeff < corrCutoff) continue;
                xic[envCol - minCol + xicStartIndex] += envelope.Abundance;
            }
            
            var smoothedXic = Smoother.Smooth(xic);
            var abundance = 0d;
            var apexScanNum = feature.MinScanNum;
            var apexIntensity = 0d;
            var boundaryIntensity = 0d;

            for (var k = 0; k < smoothedXic.Length - 1; k++)
            {
                var col = k + minCol - xicStartIndex;
                if (col < 0 || col >= NColumns - 1) continue;
                
                var centerIntensity = 0.5 * (Math.Max(0, smoothedXic[k]) + Math.Max(0, smoothedXic[k + 1]));

                if (!(centerIntensity > 0)) continue;

                var timeInterval = Run.GetElutionTime(ms1ScanNums[col + 1]) - Run.GetElutionTime(ms1ScanNums[col]);
                var abu = centerIntensity*timeInterval;
                //abuList.Add(abu);
                abundance += abu;
                if (col >= minCol && col <= maxCol && apexIntensity < abu)
                {
                    apexIntensity = abu;
                    apexScanNum = ms1ScanNums[col];
                    if (col == minCol || col == maxCol)
                    {
                        boundaryIntensity += abu;
                    }
                }
            }

            feature.SetAbundance(abundance, apexScanNum, apexIntensity, boundaryIntensity*0.5);
        }

        public double GetMs1EvidenceScore(int ms2ScanNum, double targetMass, int charge)
        {
            var tolerance = new Tolerance(Comparer.Ppm*0.5);
            if (_ms1Filter != null)
            {
                var ms2ScanNums = _ms1Filter.GetMatchingMs2ScanNums(targetMass);
                if (ms2ScanNums.Any(ms2 => ms2 == ms2ScanNum)) return 1.0d;
            }
            
            //SetTargetMass(targetMass);

            var theoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, RelativeIsotopePeakIntensityThreshold);
            var totalElutionTime = Run.GetElutionTime(Run.MaxLcScan);
            var elutionTime = Run.GetElutionTime(ms2ScanNum);
            var ms1ScanNums = Run.GetMs1ScanVector();

            var minElutionTime = elutionTime - totalElutionTime * 0.002;
            var maxElutionTime = elutionTime + totalElutionTime * 0.002;

            //var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, ms2ScanNum);
            //if (ms1ScanIndex < 0) ms1ScanIndex = ~ms1ScanIndex; // next Ms1 scan num
            var ms1ScanIndex = Run.GetPrevScanNum(ms2ScanNum, 1);

            var bcDistances = new List<double>();
            var correlations = new List<double>();
            var maxSearchScans = (int)Math.Max(ms1ScanNums.Length - ms1ScanIndex + 1, ms1ScanIndex);

            for (var i = 0; i <= maxSearchScans; i++)
            {
                for (var j = 0; j < 2; j++)
                {
                    if (i == 0 && j > 0) continue;

                    var col = (j < 1) ? ms1ScanIndex + i : ms1ScanIndex - i;
                    if (col >= ms1ScanNums.Length || col < 0) continue;
                    if (Run.GetElutionTime(ms1ScanNums[col]) > maxElutionTime || Run.GetElutionTime(ms1ScanNums[col]) < minElutionTime) continue;

                    var ms1Spec = Ms1Spectra[col];
                    var observedPeaks = ms1Spec.GetAllIsotopePeaks(targetMass, charge, theoreticalEnvelope, tolerance);

                    var bcDist = theoreticalEnvelope.GetBhattacharyyaDistance(observedPeaks);
                    var corrCoeff = theoreticalEnvelope.GetPearsonCorrelation(observedPeaks);

                    if (bcDist < 0.05)
                    {
                        return 1 - bcDist;
                    }
                    if (corrCoeff > 0.7)
                    {
                        return corrCoeff;
                    }

                    bcDistances.Add(bcDist);
                    correlations.Add(corrCoeff);
                    //enelopes.Add(obsEnv);
                }
            }

            if (bcDistances.Count < 1) return 0d;
            return correlations.Max();
        }

        public LcMsRun Run;
        public readonly List<Ms1Spectrum> Ms1Spectra;
        
        public readonly int MinSearchCharge;
        public readonly int MaxSearchCharge;
        public int MaxSearchChargeLength { get { return NRows;  } }

        public readonly int NColumns;
        public readonly int NRows;
        //public const int NRows = 30;

        public readonly MzComparerWithBinning Comparer;
        public const double RelativeIsotopePeakIntensityThreshold = 0.1d;

        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);

        private readonly int _maxThreadCount;
        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;

        private readonly double[] _distProfileAcrossCharge;
        private readonly double[] _corrProfileAcrossCharge;

        private readonly double[] _intensityAcrossCharge;

        private readonly int[,] _summedEnvelopeColRange;
        //private readonly double[] _xic;

        private LcMsPeakMatrixCell[][] _featureMatrix;
        private int[] _rows;
        private int[] _cols;
        private TheoreticalIsotopeEnvelope _theoreticalEnvelope;
        private readonly List<KeyValuePair<double, ObservedIsotopeEnvelope>> _seedEnvelopes;
        
        private void SetTargetMass(double targetMass)
        {
            /*var chargeLowerBound = Math.Ceiling(targetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Floor(targetMass / Run.MinMs1Mz);

            if (targetMass > 2000) chargeLowerBound = Math.Max(chargeLowerBound, 2);

            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);
            */
            _targetMass = targetMass;
            var chargeRange = GetDetectableMinMaxCharge(targetMass, Run.MinMs1Mz, Run.MaxMs1Mz);
            _targetMinCharge = chargeRange.Item1;
            _targetMaxCharge = chargeRange.Item2;
            //_rows = Enumerable.Range((int)rowLb, (int)(rowUb - rowLb + 1)).ToArray();
            _rows = Enumerable.Range(0, _targetMaxCharge - _targetMinCharge + 1).ToArray();
            _theoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, RelativeIsotopePeakIntensityThreshold);
        }
        
        private double _targetMass;
        private int _targetMinCharge;
        private int _targetMaxCharge;

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
            
            var observedRows        = new BitArray(NRows);
            var observedCols        = new BitArray(NColumns);
            
            var mostAbuInternalIdx = _theoreticalEnvelope.IndexOrderByRanking[0];
            var totalElutionLength = Run.GetElutionTime(Run.MaxLcScan);
            var elutionSamplingHalfLen = Math.Max(Math.Min(totalElutionLength * 0.003, 5.0), 0.5);
            var neighborHalfColumns = (int) Math.Max((elutionSamplingHalfLen/totalElutionLength)*NColumns, 5);

            var targetMassBinNum = Comparer.GetBinNumber(targetMass);
            var tolerance = new Tolerance(Comparer.Ppm*0.5);

            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;
            
            var nPeaksCutoff = NumberOfPeaksCutoff;
            var bcSeedCutoff = GetSeedBcDistThreshold();
            var corrSeedCutoff = GetSeedCorrThreshold();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            
            var options = new ParallelOptions();
            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;
            _seedEnvelopes.Clear();

            Parallel.ForEach(_rows, options, row =>
            {
                var charge = row + _targetMinCharge;

                for (var col = 0; col < NColumns; col++) _featureMatrix[row][col].Init();
               
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
                                var seed = new ObservedIsotopeEnvelope(_featureMatrix[row][col].AccurateMass, row + _targetMinCharge, ms1ScanNums[col], _featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope);

                                lock (_seedEnvelopes)
                                {
                                    _seedEnvelopes.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(bcDist, seed));
                                }
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
            if (_rows.Length < 2 || _cols.Length < 1) return clusters;

            var tempEnvelope = new double[_theoreticalEnvelope.Size];
            var tempEnvelope2 = new double[_theoreticalEnvelope.Size];
            
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            var tolerance = new Tolerance(Comparer.Ppm*0.5);

            foreach (var seed in _seedEnvelopes.OrderBy(s=>s.Key).Select(s=> s.Value))
            {
                var row = seed.Charge - _targetMinCharge;
                var col = ms1ScanNumToIndex[seed.ScanNum];
                
                if (_featureMatrix[row][col].CheckedOutFlag) continue;

                var mostAbuMz = _theoreticalEnvelope.GetIsotopeMz(seed.Charge, mostAbuInternalIndex);
                var seedLocalWin = Ms1Spectra[col].GetLocalMzWindow(mostAbuMz);
                var poissonPvalue = seedLocalWin.GetPoissonTestPvalue(_featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope.Size);
                var rankSumPvalue = seedLocalWin.GetRankSumTestPvalue(_featureMatrix[row][col].EnvelopePeaks, _theoreticalEnvelope.Size);

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

                    var minRw = (int)Math.Max(charge - _targetMinCharge - chargeNeighborGap, _rows.First());
                    var maxRw = (int)Math.Min(charge - _targetMinCharge + chargeNeighborGap, _rows.Last());
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
                            
                            if (_featureMatrix[i][j].DivergenceDist < 0.02 ||_featureMatrix[i][j].CorrelationCoeff > 0.7 || newDivergence < summedBcDist || newCorrelation > summedCorr)
                            {
                                var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass,
                                    i + _targetMinCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks,
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
                    SetCheckOutFlag(newCluster.MinCharge - _targetMinCharge, newCluster.MaxCharge - _targetMinCharge, ms1ScanNumToIndex[newCluster.MinScanNum], ms1ScanNumToIndex[newCluster.MaxScanNum], false);
                    refinedCluster = GetLcMsPeakCluster(newCluster.RepresentativeMass, newCluster.RepresentativeCharge, newCluster.MinScanNum, newCluster.MaxScanNum, true);
                }

                if (refinedCluster != null && (_scorer == null || (_scorer != null && refinedCluster.GoodEnougth && refinedCluster.Score >= _scorer.ScoreThreshold)))
                {
                    SetCheckOutFlag(_rows.First(), _rows.Last(), ms1ScanNumToIndex[refinedCluster.MinScanNum], ms1ScanNumToIndex[refinedCluster.MaxScanNum], true);
                    clusters.Add(refinedCluster);
                }
                else
                {
                    SetCheckOutFlag(newCluster.MinCharge - _targetMinCharge, newCluster.MaxCharge - _targetMinCharge, ms1ScanNumToIndex[newCluster.MinScanNum], ms1ScanNumToIndex[newCluster.MaxScanNum], true);
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

        private LcMsPeakCluster CollectLcMsPeaks(double targetMass, int minRow, int maxRow, int minCol, int maxCol, bool reCollectAllPeaks = false)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var envelopes = new List<ObservedIsotopeEnvelope>();
            var bestBcDist = 100d;
            ObservedIsotopeEnvelope bestEnvelope = null;
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            var tolerance = new Tolerance(Comparer.Ppm * 0.5);
            var massTol = tolerance.GetToleranceAsTh(targetMass);
            var nPeaksCutoff = NumberOfPeaksCutoff;

            var bcCutoff = GetSeedBcDistThreshold();
            var corrCutoff = GetSeedCorrThreshold();

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (reCollectAllPeaks) _featureMatrix[i][j].Init();

                    if (reCollectAllPeaks || !_featureMatrix[i][j].Exist || Math.Abs(_featureMatrix[i][j].AccurateMass - targetMass) > massTol)
                    {
                        var peaks = Ms1Spectra[j].GetAllIsotopePeaks(targetMass, i + _targetMinCharge, _theoreticalEnvelope, tolerance);

                        if (peaks.Count(p => p != null) > 0)
                        {
                            _featureMatrix[i][j].DivergenceDist = _theoreticalEnvelope.GetBhattacharyyaDistance(peaks); ;
                            _featureMatrix[i][j].AccurateMass = targetMass;
                            _featureMatrix[i][j].CorrelationCoeff = _theoreticalEnvelope.GetPearsonCorrelation(peaks); ;
                            Array.Copy(peaks, _featureMatrix[i][j].EnvelopePeaks, peaks.Length);    
                        }
                    }

                    if (!_featureMatrix[i][j].Exist) continue;
                    if (_featureMatrix[i][j].CountActivePeaks < nPeaksCutoff) continue;
                    if (_featureMatrix[i][j].DivergenceDist > bcCutoff && _featureMatrix[i][j].CorrelationCoeff < corrCutoff) continue; // exclude outliers
                    var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + _targetMinCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, _theoreticalEnvelope);
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
            cluster.AddEnvelopes(minRow + _targetMinCharge, maxRow + _targetMinCharge, ms1ScanNums[minCol], ms1ScanNums[maxCol], envelopes);
            
            return cluster;
        }

        public LcMsPeakCluster CollectLcMsPeaksWithNoise(double targetMass, int targetCharge, int targetMinScanNum, int targetMaxScanNum, int targetMinCharge, int targetMaxCharge)
        {
            SetTargetMass(targetMass);

            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            if (Run.GetMsLevel(targetMinScanNum) > 1) targetMinScanNum = Run.GetPrevScanNum(targetMinScanNum, 1);
            if (Run.GetMsLevel(targetMaxScanNum) > 1) targetMaxScanNum = Run.GetNextScanNum(targetMaxScanNum, 1);

            var minCol = ms1ScanNumToIndex[targetMinScanNum];
            var maxCol = ms1ScanNumToIndex[targetMaxScanNum];

            var minRow = targetMinCharge - _targetMinCharge;
            var maxRow = targetMaxCharge - _targetMinCharge;
            var abundance = 0d;

            for (var j = minCol; j <= maxCol; j++)
            {
                //var mostAbuMz = _theoreticalEnvelope.GetIsotopeMz(i + MinScanCharge, mostAbuInternalIndex);
                //var seedLocalWin = Ms1Spectra[j].GetLocalMzWindow(mostAbuMz);
                //_xic[j] += seedLocalWin.MedianIntensity;
                //tempAbundance += Ms1Spectra[j].MedianIntensity;
                abundance += Ms1Spectra[j].MedianIntensity;
            }

            /*
            Array.Clear(_xic, 0, _xic.Length);
            var tempAbundance = 0d;
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist) continue;
                    _xic[j] += _featureMatrix[i][j].Intensity;
                    tempAbundance += _featureMatrix[i][j].Intensity;
                }
            }

            if (!(tempAbundance > 0))
            {
                for (var i = minRow; i <= maxRow; i++)
                {
                    for (var j = minCol; j <= maxCol; j++)
                    {
                        var mostAbuMz = _theoreticalEnvelope.GetIsotopeMz(i + MinScanCharge, mostAbuInternalIndex);
                        var seedLocalWin = Ms1Spectra[j].GetLocalMzWindow(mostAbuMz);
                        _xic[j] += seedLocalWin.MedianIntensity;
                        tempAbundance += Ms1Spectra[j].MedianIntensity;
                    }
                }
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
            
            if (!(abundance > 0)) abundance = tempAbundance;
            */

            var repScanNum = ms1ScanNums[(int)((minCol+maxCol)*0.5)];
            var repMz = 0;

            var cluster = new LcMsPeakCluster(Run, _theoreticalEnvelope, targetMass, targetCharge, repMz, repScanNum, abundance);
            cluster.AddEnvelopes(minRow + _targetMinCharge, maxRow + _targetMinCharge, ms1ScanNums[minCol], ms1ScanNums[maxCol]);
            cluster.Score = -999d;

            return cluster;
        }


        private double GetBcDistThreshold()
        {
            if (_targetMass > 35000) return 0.12;
            if (_targetMass > 25000) return 0.1;
            if (_targetMass > 15000) return 0.08;
            if (_targetMass > 5000) return 0.06;
            
            return 0.03;
        }

        private double GetCorrThreshold()
        {
            if (_targetMass > 35000) return 0.3;
            if (_targetMass > 25000) return 0.4;
            if (_targetMass > 15000) return 0.5;
            if (_targetMass > 5000) return 0.6;
            return 0.7;
        }
        
        private double GetSeedBcDistThreshold()
        {
            if (_targetMass > 45000) return 0.3;
            if (_targetMass > 35000) return 0.25;
            if (_targetMass > 25000) return 0.2;
            if (_targetMass > 15000) return 0.15;
            if (_targetMass > 10000) return 0.1;
            if (_targetMass > 3000) return 0.07;;
            return 0.05;
        }

        private double GetSeedCorrThreshold()
        {
            if (_targetMass > 45000) return 0.2;
            if (_targetMass > 35000) return 0.3;
            if (_targetMass > 25000) return 0.4;
            if (_targetMass > 15000) return 0.5;
            if (_targetMass > 10000) return 0.6;
            if (_targetMass > 3000) return 0.7;
            return 0.8;
        }

        private Tuple<int, int, int, int> DetermineFeatureRange(double targetMass, int seedRow, int seedMinCol, int seedMaxCol)
        {
            Array.Clear(_distProfileAcrossCharge, 0, _distProfileAcrossCharge.Length);
            Array.Clear(_corrProfileAcrossCharge, 0, _corrProfileAcrossCharge.Length);
            Array.Clear(_intensityAcrossCharge, 0, _intensityAcrossCharge.Length);

            if (_rows.Length < 2 || _cols.Length < 1) return null;

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

                _intensityAcrossCharge[i] = summedEnvelope.Sum();
            });
            
            foreach (var i in _rows)
            {
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

                var charge = i + _targetMinCharge;
                var chargeIdx = (charge % 2 == 0) ? 0 : 1;
                if (bestCharge[chargeIdx] == 0 || _distProfileAcrossCharge[i] < bestChargeDist[chargeIdx])
                {
                    bestChargeDist[chargeIdx] = _distProfileAcrossCharge[i];
                    bestCharge[chargeIdx] = charge;
                }
            }
            
            if (bestCharge[0] == 0 && bestCharge[1] == 0) return null;

            // only one charge? there should be another....force to cover neighboring charge states
            if (bestCharge[0] == 0 || bestCharge[1] == 0)
            {
                minRow = Math.Max(minRow - 1, _rows.First());
                maxRow = Math.Min(maxRow + 1, _rows.Last());
            }

            if (minRow == maxRow) return null;

            for (var i = minRow; i <= maxRow; i++)
            {
                var charge = i + _targetMinCharge;
                var chargeIdx = (charge % 2 == 0) ? 0 : 1;
                if (bestCharge[chargeIdx] == 0)
                {
                    bestCharge[chargeIdx] = charge;
                    break;
                }
            }

            // refine scan boundary by considering bell shaped Xic for the best charge state
            //Array.Clear(_xic, 0, _xic.Length);

            var xic = new double[maxCol - minCol + 1 + 18];
            const int xicStartIndex = 9;
            var xicEndIndex = xic.Length - 9 - 1;

            var nIntensity = 0;
            foreach (var c in bestCharge)
            {
                var i = c - _targetMinCharge;
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist ||
                        (_featureMatrix[i][j].DivergenceDist > GetSeedBcDistThreshold() && _featureMatrix[i][j].CorrelationCoeff < GetSeedCorrThreshold())) continue;
                    //_xic[j] += _featureMatrix[i][j].Intensity;

                    xic[j - minCol + xicStartIndex] += _featureMatrix[i][j].Intensity;
                    nIntensity++;
                }
            }

            if (nIntensity < 2)
            {
                return null;
            }

            var smoothedXic = Smoother.Smooth(xic);
            var apexCol = -1;
            var apexIntensity = 0d;
            
            //for (var j = minCol; j <= maxCol; j++)
            for (var j = xicStartIndex; j <= xicEndIndex; j++)
            {
                if (smoothedXic[j] > apexIntensity)
                {
                    apexIntensity = smoothedXic[j];
                    apexCol = j - xicStartIndex + minCol;
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
                if (smoothedXic[j - minCol + xicStartIndex] < oneSigIntensity) break;
            }

            var oneSigPeriod = Run.GetElutionTime(ms1ScanNums[apexCol]) -
                               Run.GetElutionTime(ms1ScanNums[elutionStartColByOneSigma]);

            // find 2-sigma position at elution time prior to Apex (going backward)
            /*var elutionStartColByTwoSigma = elutionStartColByOneSigma;
            for (var j = elutionStartColByOneSigma - 1; j >= 0; j--)
            {
                elutionStartColByTwoSigma = j;
                if (Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]) > 2*oneSigPeriod) break;
                var needBreak = true;
                for (var k = j - 1; k >= j - 3; k--)
                {
                    if (k < 0) break;
                    if (smoothedXic[j] * 1.3 < smoothedXic[k])
                    {
                        needBreak = false;
                        break;
                    }
                }
                if (needBreak) break;
            }*/

            // extends for long tail
            //var elutionStartCol = elutionStartColByTwoSigma;
            var elutionStartCol = elutionStartColByOneSigma;
            //for (var j = elutionStartColByTwoSigma - 1; j >= minCol; j--)
            for (var j = elutionStartColByOneSigma - 1; j >= minCol; j--)
            {
                //if (smoothedXic[j - minCol + xicStartIndex] < threeSigIntensity) break;
                if (smoothedXic[j - minCol + xicStartIndex] <= double.Epsilon && xic[j - minCol + xicStartIndex] <= double.Epsilon) break;
                var elutionLen = Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]);
                //if (elutionLen > 3 && elutionLen > 4 * oneSigPeriod) break;

                var needBreak = true;
                for (var k = j - 1; k >= j - 3; k--)
                {
                    if (k < 0) break;
                    if (smoothedXic[k - minCol + xicStartIndex] < smoothedXic[j - minCol + xicStartIndex] * 1.3)
                    {
                        needBreak = false;
                        break;
                    }
                }
                if (needBreak) break;
                    
                elutionStartCol = j;
            }

            ////////////////// Second Half
            var elutionEndColByOneSigma = apexCol;
            for (var j = apexCol + 1; j < NColumns; j++)
            {
                elutionEndColByOneSigma = j;
                if (smoothedXic[j - minCol + xicStartIndex] < oneSigIntensity) break;
            }
            // find two sigma point at elution time axis
            /*
            oneSigPeriod = Run.GetElutionTime(ms1ScanNums[elutionEndColByOneSigma]) - Run.GetElutionTime(ms1ScanNums[apexCol]);
            var elutionEndColByTwoSigma = elutionEndColByOneSigma;
            for (var j = elutionEndColByOneSigma + 1; j < NColumns; j++)
            {
                elutionEndColByTwoSigma = j;
                if (Run.GetElutionTime(ms1ScanNums[j]) -  Run.GetElutionTime(ms1ScanNums[apexCol]) > 2 * oneSigPeriod) break;
                
                var needBreak = true;
                for (var k = j + 1; k <= j + 3; k++)
                {
                    if (k >= NColumns) break;
                    if (smoothedXic[j] * 1.3 > smoothedXic[k])
                    {
                        needBreak = false;
                        break;
                    }
                }
                if (needBreak)
                {
                    break;
                }
            }*/

            // possible extension?
            //var elutionEndCol = elutionEndColByTwoSigma;
            var elutionEndCol = elutionEndColByOneSigma;
            //for (var j = elutionEndColByTwoSigma + 1; j <= maxCol; j++)
            for (var j = elutionEndColByOneSigma + 1; j <= maxCol; j++)
            {
                //if (smoothedXic[j - minCol + xicStartIndex] < threeSigIntensity) break;
                if (smoothedXic[j - minCol + xicStartIndex] <= double.Epsilon && xic[j - minCol + xicStartIndex] <= double.Epsilon) break;

                var elutionLen = Run.GetElutionTime(ms1ScanNums[j]) - Run.GetElutionTime(ms1ScanNums[apexCol]);
                if (elutionLen > 3 && elutionLen > 8 * oneSigPeriod) break;

                
                    var needBreak = true;
                    for (var k = j + 1; k <= j + 3; k++)
                    {
                        if (k >= NColumns) break;
                        if (smoothedXic[j - minCol + xicStartIndex] * 1.3 > smoothedXic[k - minCol + xicStartIndex])
                        {
                            needBreak = false;
                            break;
                        }
                    }
                    if (needBreak)
                    {
                        break;
                    }

                elutionEndCol = j;
            }
            
            return new Tuple<int, int, int, int>(minRow, maxRow, elutionStartCol, elutionEndCol);
        }

        private double[] GetSummedEnvelopeAtCharge(double targetMass, int row, int minCol, int maxCol)
        {
            const int maxScanSkips = 2;
            const double goodEnoughBcDistance = 0.07;
            const double goodEnoughCorrCoeff = 0.7;
            var summedEnvelope = new double[_theoreticalEnvelope.Size];
            var seedCol = -1;
            var seedDist = 10.0d;
            
            var newMinCol = minCol;
            var newMaxCol = maxCol;
            //var tolerance = new Tolerance(5);
            var tolerance = new Tolerance(Comparer.Ppm * 0.5);
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

            if (envelope.NumberOfPeaks > nPeaks * 0.7) return true;

            //var tolerance = new Tolerance(5);
            var tolerance = new Tolerance(Comparer.Ppm * 0.5);
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

        public Tuple<int, int> GetDetectableMinMaxCharge(double mass, double minMs1Mz, double maxMs1Mz)
        {
            if (mass < 1000) return new Tuple<int, int>(MinSearchCharge, Math.Min(MinSearchCharge + MaxSearchChargeLength - 1, 5));
            if (mass < 2500) return new Tuple<int, int>(MinSearchCharge, Math.Min(MinSearchCharge + MaxSearchChargeLength - 1, 10));
            if (mass > 40000) return new Tuple<int, int>(MaxSearchCharge - MaxSearchChargeLength + 1, MaxSearchCharge);
            
            const double c1 = 8.52456367162247e-09;
            const double c2 = 0.000389353356727307;
            const double c3 = -1.50236827038509;
            const double p1 = -1.29109336042109e-08;
            const double p2 = 0.00180628243605134;
            const double p3 = 9.35391578169730;

            var chargeMin = (int)Math.Round((c1 * mass) * (c1 * mass) + c2 * mass + c3);
            var chargeMax = (int) Math.Round((p1 * mass) * (p1 * mass) + p2 * mass + p3);

            var chargeUb = (int)Math.Min(Math.Floor(mass / minMs1Mz), MaxSearchCharge);
            var chargeLb = (int) Math.Max(Math.Ceiling(mass/maxMs1Mz), MinSearchCharge);

            var minCharge = (int) Math.Max(chargeLb, chargeMin);
            var maxCharge = (int) Math.Min(chargeUb, chargeMax);

            if (maxCharge - minCharge + 1 > MaxSearchChargeLength)
            {
                maxCharge = minCharge + MaxSearchChargeLength - 1;
            }

            return new Tuple<int, int>(minCharge, maxCharge);
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
