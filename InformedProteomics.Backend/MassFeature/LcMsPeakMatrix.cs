using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security.Principal;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MSFileReaderLib;
using Constants = InformedProteomics.Backend.Data.Biology.Constants;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakMatrix
    {
        public LcMsPeakMatrix(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0, LcMsFeatureLikelihood scorer = null)
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

            _distProfileAcrossCharge = new double[NRows];
            _summedEnvelopeColRange = new int[NRows, 2];
            
            _xic = new double[NColumns];
            _featureMatrix = null;
            Comparer = new MzComparerWithBinning(28);
            _scorer = scorer;
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
        /*
        public void EvaludateSummedDecoyEnvelope(LcMsPeakCluster feature)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var minCol = ms1ScanNumToIndex[feature.MinScanNum];
            var maxCol = ms1ScanNumToIndex[feature.MaxScanNum];

            var chargeLowerBound = Math.Floor(feature.Mass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Ceiling(feature.Mass / Run.MinMs1Mz);
            var minRow = (int)Math.Max(0, chargeLowerBound - MinScanCharge);
            var maxRow = (int)Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);

            int n = _ms1PeakList.Count;
            var rnd = new Random();

            var summedPeakCount = new int[NRows];
            var summedPeakIntensity = new double[NRows];
            var summedMedianIntensity = new double[NRows];
            var nRows = maxRow - minRow + 1;

            feature.EnvelopeDistanceScoreAcrossCharge = new double[nRows];
            feature.EnvelopeIntensityScoreAcrossCharge = new double[nRows];
            feature.AbundanceDistributionAcrossCharge = new double[nRows];

            //AbundanceDistributionAcrossCharge = new double[nRows];

            feature.BestCorrelationScoreAcrossCharge = new double[nRows];
            feature.BestDistanceScoreAcrossCharge = new double[nRows];
            feature.BestIntensityScoreAcrossCharge = new double[nRows];
            
            EnvelopeDistanceScoreAcrossCharge = new double[nRows];
            feature.EnvelopeCorrelationScoreAcrossCharge = new double[nRows];
            //EnvelopeIntensityScoreAcrossCharge = new double[nRows];

            var intensities = new double[_theoreticalEnvelope.Size];

            var comparer = new MzComparerWithBinning(28);
            var numMzBins = comparer.GetBinNumber(Run.MaxMs1Mz) - comparer.GetBinNumber(Run.MinMs1Mz) + 1;

            for (var i = minRow; i <= maxRow; i++)
            {
                Array.Clear(intensities, 0, intensities.Length);
                Array.Clear(summedPeakCount, 0, summedPeakCount.Length);
                Array.Clear(summedPeakIntensity, 0, summedPeakIntensity.Length);
                Array.Clear(summedMedianIntensity, 0, summedMedianIntensity.Length);

                for (var j = minCol; j <= maxCol; j++)
                {
                    for (var k = 0; k < _theoreticalEnvelope.Size; k++)
                    {
                        var r = rnd.Next(0, numMzBins);
                        if (r < Ms1Spectra[j].Peaks.Length)
                        {
                            r = rnd.Next(0, n);
                            var randomPeak = _ms1PeakList[r];
                            intensities[k] += randomPeak.Intensity;

                            if (k == _theoreticalEnvelope.IndexOrderByRanking[0])
                            {
                                summedPeakCount[i]++;
                                summedPeakIntensity[i] += randomPeak.Intensity;
                                summedMedianIntensity[i] += Ms1Spectra[randomPeak.Ms1SpecIndex].GetLocalMedianIntensity(randomPeak, feature.Mass);
                            }
                        }

                    }
                }

                feature.EnvelopeDistanceScoreAcrossCharge[i - minRow] = _theoreticalEnvelope.GetBhattacharyyaDistance(intensities);
                if (summedMedianIntensity[i] > 0) feature.EnvelopeIntensityScoreAcrossCharge[i - minRow] = summedPeakIntensity[i] / summedMedianIntensity[i];
                feature.AbundanceDistributionAcrossCharge[i - minRow] = intensities.Sum();
            }

            var s = feature.AbundanceDistributionAcrossCharge.Sum();
            foreach (var i in _rows)
                feature.AbundanceDistributionAcrossCharge[i - minRow] = feature.AbundanceDistributionAcrossCharge[i - minRow] / s;

            feature.SetChargeRange(minRow + MinScanCharge, maxRow + MinScanCharge);
        }*/

        public LcMsPeakCluster GetLcMsPeakCluster(double targetMass, int targetCharge, int targetMinScanNum, int targetMaxScanNum, int seedMs1ScanNum = -1, bool featureMatrixCreated = false)
        {
            if (!featureMatrixCreated) BuildFeatureMatrix(targetMass); // should be called first

            // 1) set initial elution period for initial sampling
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            if (Run.GetMsLevel(targetMinScanNum) > 1) targetMinScanNum = Run.GetPrevScanNum(targetMinScanNum, 1);
            if (Run.GetMsLevel(targetMaxScanNum) > 1) targetMaxScanNum = Run.GetNextScanNum(targetMaxScanNum, 1);

            var seedCol = -1;
            var row = Math.Max(Math.Min(targetCharge - MinScanCharge, MaxScanCharge - MinScanCharge), 0);

            if (seedMs1ScanNum < 0)
            {
                var targetScanNum = (int)((targetMinScanNum + targetMaxScanNum) * 0.5d);
                seedCol = Array.BinarySearch(ms1ScanNums, targetScanNum);
                if (seedCol < 0) seedCol = ~seedCol;
            }
            else
            {
                seedCol = ms1ScanNumToIndex[seedMs1ScanNum];
            }

            var elutionSamplingPeriod = Math.Max(Math.Min(Run.GetElutionTime(Run.MaxLcScan) * 0.003, 5.0), 0.5);
            var tempWindow = GetElutionWindow(seedCol, elutionSamplingPeriod);
            var minCol = Math.Min(tempWindow.Item1, ms1ScanNumToIndex[targetMinScanNum]);
            var maxCol = Math.Max(tempWindow.Item2, ms1ScanNumToIndex[targetMaxScanNum]);

            // 2) determine charge state range and elution period
            
            var range = DetermineFeatureRange(row, minCol, maxCol, seedCol);
            var minRow = range.Item1;
            var maxRow = range.Item2;
            minCol = range.Item3;
            maxCol = range.Item4;

            // 3) collect envelopes
            var feature = CollectLcMsPeaks(minRow, maxRow, minCol, maxCol);

            if (feature == null) return null;

            feature.UpdateScore(Ms1Spectra);
            if (_scorer != null) feature.Score = _scorer.GetScore(feature);

            // 4) determine abundance
            Array.Clear(_xic, 0, _xic.Length);
            foreach (var envelope in feature.Envelopes)
            {
                // skip charge states having no good summed envelope
                if (feature.EnvelopeDistanceScoreAcrossCharge[envelope.Charge - feature.MinCharge] > 0.1) continue;
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
        /*
        public double[] GetEntireXic(LcMsPeakCluster feature)
        {
            Array.Clear(_xic, 0, _xic.Length);
            var row = feature.RepresentativeCharge - MinScanCharge;
            foreach (var col in _cols)
            {
                _xic[col] = _featureMatrix[row][col].Intensity;
            }
            var smoothedXic = Smoother.Smooth(_xic);
            return smoothedXic;
        }*/


        public LcMsRun Run;
        public readonly List<Ms1Spectrum> Ms1Spectra;
        public readonly int MinScanCharge;
        public readonly int MaxScanCharge;
        public readonly Tolerance MzTolerance;

        public readonly int NColumns;
        public readonly int NRows;
        public readonly MzComparerWithBinning Comparer;
        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);

        private readonly int _maxThreadCount;
        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;

        private readonly double[] _distProfileAcrossCharge;
        private readonly int[,] _summedEnvelopeColRange;
        private readonly double[] _xic;

        private LcMsPeakMatrixCell[][] _featureMatrix;
        private int[] _rows;
        private int[] _cols;
        
        private TheoreticalIsotopeEnvelope _theoreticalEnvelope;
        
        
        private void SetTargetMass(double targetMass)
        {
            var chargeLowerBound = Math.Floor(targetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Ceiling(targetMass / Run.MinMs1Mz);
            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);
            _targetMass = targetMass;
            _rows = Enumerable.Range((int)rowLb, (int)(rowUb - rowLb + 1)).ToArray();
            _theoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, 0.1d);
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

        private void BuildFeatureMatrix(double targetMass)
        {
            InitFeatureMatrix();

            SetTargetMass(targetMass);

            //var queryMassBinNum = Comparer.GetBinNumber(targetMass);

            var options             = new ParallelOptions();
            var observedRows        = new bool[NRows];
            var observedCols        = new bool[NColumns];

            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            //var mostAbuIsotopeInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];
            //var mostAbuIsotopeIndex = _theoreticalEnvelope.GetMostAbundantIsotope().Index;
            
            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;

            const double signalToNoiseRatioCutoff = 1.4826;
            var nPeaksCutoff = (targetMass > 2000) ? 3 : 2;

            Parallel.ForEach(_rows, options, row =>
            {
                for (var col = 0; col < NColumns; col++) _featureMatrix[row][col].Init();

                var charge = row + MinScanCharge;
                for (var k = 0; k < _theoreticalEnvelope.Size; k++)
                {
                    var i = _theoreticalEnvelope.IndexOrderByRanking[k];
                    var isotopeIndex = _theoreticalEnvelope.Isotopes[i].Index;

                    
                    var isotopeMz = Ion.GetIsotopeMz(targetMass, charge, isotopeIndex);
                    var mzTol = MzTolerance.GetToleranceAsTh(isotopeMz);
                    var isotopeMzLb = isotopeMz - mzTol;
                    var isotopeMzUb = isotopeMz + mzTol;

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
                            if (!(_featureMatrix[row][col].AccurateMass > 0)) _featureMatrix[row][col].AccurateMass = targetMass;
                          
                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (_featureMatrix[row][col].EnvelopePeaks[i] != null)
                            {
                                var tmpPeak = _featureMatrix[row][col].EnvelopePeaks[i];
                                var bc1 = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                                _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                var bc2 = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                                if (bc1 < bc2) _featureMatrix[row][col].EnvelopePeaks[i] = tmpPeak;
                            }
                            else
                            {
                                _featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                            }
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
                        _featureMatrix[row][col].CorrelationCoeff = _theoreticalEnvelope.GetPearsonCorrelation(_featureMatrix[row][col].EnvelopePeaks);
                        _featureMatrix[row][col].DivergenceDist = _theoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[row][col].EnvelopePeaks);
                        if (!observedRows[row]) observedRows[row] = true;
                        if (!observedCols[col]) observedCols[col] = true;
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

        private IEnumerable<ObservedIsotopeEnvelope> GetSeedCells(double targetMass)
        {
            var seedList = new List<KeyValuePair<double, ObservedIsotopeEnvelope>>();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];

            var nPeaksCutoff = 2;
            if (targetMass > 25000) nPeaksCutoff = 5;
            else if (targetMass > 8000) nPeaksCutoff = 4;
            else if (targetMass > 2000) nPeaksCutoff = 3;
            else nPeaksCutoff = 2;

            foreach (var i in _rows)
            {
                foreach (var j in _cols)
                {
                    if (!_featureMatrix[i][j].Exist) continue;
                    var mostAbuPeak = _featureMatrix[i][j].EnvelopePeaks[mostAbuInternalIndex];
                    if (mostAbuPeak == null) continue;

                    var signalToNoiseRatio = mostAbuPeak.Intensity / Ms1Spectra[j].MedianIntensity;
                    if (signalToNoiseRatio < 3) continue;

                    if (_featureMatrix[i][j].EnvelopePeaks.Count(p => p != null && p.Active) < nPeaksCutoff) continue;

                    var bcDist = _featureMatrix[i][j].DivergenceDist;
                    var corr = _featureMatrix[i][j].CorrelationCoeff;

                    if (bcDist > 0.2 || corr < 0.3) continue;

                    var seed = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, _theoreticalEnvelope);
                    seedList.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(bcDist, seed));
                }
            }
            return seedList.OrderBy(x => x.Key).Select(x => x.Value);
        }

        private IList<LcMsPeakCluster> GetLcMs1PeakClusters(int binNumber)
        {
            var targetMass = Comparer.GetMzAverage(binNumber);
            BuildFeatureMatrix(targetMass); // should be called first

            var clusters = new List<LcMsPeakCluster>();
            var tempEnvelope = new double[_theoreticalEnvelope.Size];
            var tempEnvelope2 = new double[_theoreticalEnvelope.Size];
            
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            foreach (var seed in GetSeedCells(targetMass))
            {
                var row = seed.Charge - MinScanCharge;
                var col = ms1ScanNumToIndex[seed.ScanNum];

                if (_featureMatrix[row][col].CheckedOutFlag) continue;
                var seedMass = _featureMatrix[row][col].AccurateMass;
                var massTol = MzTolerance.GetToleranceAsTh(seedMass);
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
                    var chargeNeighborGap = GetScatteredChargeLength(charge);
                    var minRw = (int) Math.Max(charge - MinScanCharge - chargeNeighborGap, _rows.First());
                    var maxRw = (int) Math.Min(charge - MinScanCharge + chargeNeighborGap, _rows.Last());

                    var currRow = charge - MinScanCharge;
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
                            
                            if (_featureMatrix[i][j].DivergenceDist < 0.04 &&
                                _featureMatrix[i][j].CorrelationCoeff > 0.8 || newDivergence < summedBcDist ||
                                newCorrelation > summedCorr)
                            {
                                var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass,
                                    i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks,
                                    _theoreticalEnvelope);

                                neighbors.Enqueue(envelope);
                                newCluster.AddObservedEnvelope(envelope);
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

                if (summedCorr > 0.3 || summedBcDist < 0.2)
                {
                    var refinedCluster = GetLcMsPeakCluster(newCluster.RepresentativeMass, newCluster.RepresentativeCharge, newCluster.MinScanNum, newCluster.MaxScanNum, newCluster.RepresentativeScanNum, true);
                    if (refinedCluster != null && refinedCluster.Score > -10)
                    {
                        // re-update check-out map
                        for (var i = newCluster.MinCharge - MinScanCharge;
                            i <= newCluster.MaxCharge - MinScanCharge;
                            i++)
                            for (var j = ms1ScanNumToIndex[newCluster.MinScanNum];
                                j <= ms1ScanNumToIndex[newCluster.MaxScanNum];
                                j++) _featureMatrix[i][j].CheckedOutFlag = false;

                        for (var i = refinedCluster.MinCharge - MinScanCharge;
                            i <= refinedCluster.MaxCharge - MinScanCharge;
                            i++)
                            for (var j = ms1ScanNumToIndex[refinedCluster.MinScanNum];
                                j <= ms1ScanNumToIndex[refinedCluster.MaxScanNum];
                                j++) _featureMatrix[i][j].CheckedOutFlag = true;


                        refinedCluster.tempInitialCorr = summedCorr;
                        refinedCluster.tempInitialDist = summedBcDist;
                        clusters.Add(refinedCluster);
                    }
                    else
                    {
                        for (var i = newCluster.MinCharge - MinScanCharge; i <= newCluster.MaxCharge - MinScanCharge; i++)
                            for (var j = ms1ScanNumToIndex[newCluster.MinScanNum]; j <= ms1ScanNumToIndex[newCluster.MaxScanNum]; j++) _featureMatrix[i][j].CheckedOutFlag = true;                                            
                    }
                }
                else
                {
                    for (var i = newCluster.MinCharge - MinScanCharge; i <= newCluster.MaxCharge - MinScanCharge; i++)
                        for (var j = ms1ScanNumToIndex[newCluster.MinScanNum]; j <= ms1ScanNumToIndex[newCluster.MaxScanNum]; j++) _featureMatrix[i][j].CheckedOutFlag = true;                    
                }
            }
            return clusters;            
        }

        private int GetScatteredChargeLength(int charge)
        {
            /*
            var scatteredLength = 2;
            if (charge > 40) scatteredLength = 12;
            else if (charge > 30) scatteredLength = 9;
            else if (charge > 20) scatteredLength = 6;
            else if (charge > 10) scatteredLength = 3;
            return scatteredLength;
            */
            return 3;
        }

        private LcMsPeakCluster CollectLcMsPeaks(int minRow, int maxRow, int minCol, int maxCol)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var envelopes = new List<ObservedIsotopeEnvelope>();
            var bestBcDist = 100d;
            ObservedIsotopeEnvelope bestEnvelope = null;
            var mostAbuInternalIndex = _theoreticalEnvelope.IndexOrderByRanking[0];

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_featureMatrix[i][j].Exist) continue;

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
            foreach(var e in envelopes) cluster.AddObservedEnvelope(e);
            
            return cluster;
        }

        private Tuple<int, int, int, int> DetermineFeatureRange(int seedRow, int seedMinCol, int seedMaxCol, int seedCol)
        {
            var bestBcDist = 10.0d;
            var bestBcDistRow = 0;
            Array.Clear(_distProfileAcrossCharge, 0, _distProfileAcrossCharge.Length);
            
            // for each carge state, find and evaluate optimial summed envelope
            foreach (var i in _rows)
            {
                var c1 = 0;
                var c2 = 0 ;
                var summedBcDist = GetSummedEnvelopeAtCharge(i, seedMinCol, seedMaxCol, out c1, out c2);
                _distProfileAcrossCharge[i] = summedBcDist;

                if (c1 <= seedCol && seedCol <= c2 && _distProfileAcrossCharge[i] < bestBcDist)
                {
                    bestBcDistRow = i;
                    bestBcDist = _distProfileAcrossCharge[i];
                }

                _summedEnvelopeColRange[i,0] = c1;
                _summedEnvelopeColRange[i,1] = c2;
            }

            //var bcCutoff = OtsuThreshold.GetThreshold(_distProfileAcrossCharge, 0, 0.3, 0.0025, 1e-9, 0.3);
            //bcCutoff = Math.Min(Math.Max(bcCutoff, 0.05), 0.1);
            const double bcCutoff = 0.1;

            var minCol = NColumns - 1;
            var maxCol = 0;

            var minColBestCharge = _summedEnvelopeColRange[bestBcDistRow,0];
            var maxColBestCharge = _summedEnvelopeColRange[bestBcDistRow,1];

            // determine temporary charge and scan boundaries
            var minRow = _rows.Last();
            var maxRow = _rows.First();
            foreach (var i in _rows)
            {
                if (_distProfileAcrossCharge[i] > bcCutoff) continue;

                if (_summedEnvelopeColRange[i,1] < minColBestCharge) continue;
                if (_summedEnvelopeColRange[i,0] > maxColBestCharge) continue;

                if (i < minRow) minRow = i;
                if (i > maxRow) maxRow = i;
                
                if (_summedEnvelopeColRange[i,0] < minCol) minCol = _summedEnvelopeColRange[i,0];
                if (_summedEnvelopeColRange[i,1] > maxCol) maxCol = _summedEnvelopeColRange[i,1];
            }

            // refine scan boundary by considering bell shaped Xic for the best charge state
            Array.Clear(_xic, 0, _xic.Length);
            for (var j = minCol; j <= maxCol; j++) _xic[j] = _featureMatrix[bestBcDistRow][j].Intensity;

            var smoothedXic = Smoother.Smooth(_xic);
            var apexCol = 0;
            var apexIntensity = 0d;
            for (var j = minColBestCharge; j <= maxColBestCharge; j++)
            {
                if (smoothedXic[j] > apexIntensity)
                {
                    apexIntensity = smoothedXic[j];
                    apexCol = j;
                }
            }
            
            ////////////////// 1) First Half
            var ms1ScanNums = Run.GetMs1ScanVector();
            var oneSigIntensity = apexIntensity * Math.Exp(-1);
            var threeSigIntensity = apexIntensity * Math.Exp(-4.5);
            
            // estimate sigma for Gaussian shaped elution profile for the first half
            var elutionStartColByOneSigma = 0;
            for (var j = apexCol - 1; j >= 0; j--)
            {
                if (smoothedXic[j] < oneSigIntensity)
                {
                    elutionStartColByOneSigma = j;
                    break;
                }
            }

            var oneSigPeriod = Run.GetElutionTime(ms1ScanNums[apexCol]) -
                               Run.GetElutionTime(ms1ScanNums[elutionStartColByOneSigma]);

            // find 2-sigma position at elution time prior to Apex (going backward)
            var elutionStartColByTwoSigma = elutionStartColByOneSigma;
            for (var j = elutionStartColByOneSigma - 1; j >= 0; j--)
            {
                if (Run.GetElutionTime(ms1ScanNums[j]) < Run.GetElutionTime(ms1ScanNums[apexCol]) - 2*oneSigPeriod)
                    elutionStartColByTwoSigma = j;
                else break;
            }

            // extends for long tail
            var elutionStartCol = elutionStartColByTwoSigma;
            for (var j = elutionStartColByTwoSigma - 1; j >= minCol; j--)
            {
                if (smoothedXic[j] < threeSigIntensity) break;
                if (Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]) > oneSigPeriod*3) break;

                var needBreak = true;
                for (var k = j; k >= j - 3; k--)
                {
                    if (k < minCol) break;
                    if (smoothedXic[k] < smoothedXic[j + 1]) needBreak = false;
                }
                if (needBreak) break;
                    
                elutionStartCol = j;
            }

            ////////////////// Second Half
            var elutionEndColByOneSigma = 0;
            for (var j = apexCol + 1; j < NColumns; j++)
            {
                if (smoothedXic[j] < oneSigIntensity)
                {
                    elutionEndColByOneSigma = j;
                    break;
                }
            }
            // find two sigma point at elution time axis
            oneSigPeriod = Run.GetElutionTime(ms1ScanNums[elutionEndColByOneSigma]) - Run.GetElutionTime(ms1ScanNums[apexCol]);
            var elutionEndColByTwoSigma = elutionEndColByOneSigma;
            for (var j = elutionEndColByOneSigma + 1; j < NColumns; j++)
            {
                if (Run.GetElutionTime(ms1ScanNums[j]) < Run.GetElutionTime(ms1ScanNums[apexCol]) + 2 * oneSigPeriod)
                    elutionEndColByTwoSigma = j;
                else break;
            }

            // possible extension?
            var elutionEndCol = elutionEndColByTwoSigma;
            for (var j = elutionEndColByTwoSigma + 1; j <= maxCol; j++)
            {
                if (smoothedXic[j] < threeSigIntensity) break;
                if (Run.GetElutionTime(ms1ScanNums[apexCol]) - Run.GetElutionTime(ms1ScanNums[j]) > oneSigPeriod * 4) break;

                var needBreak = true;
                for (var k = j; k <= j + 3; k++)
                {
                    if (k > maxCol) break;
                    if (smoothedXic[k] < smoothedXic[j - 1]) needBreak = false;
                }
                if (needBreak) break;

                elutionEndCol = j;
            }
            
            return new Tuple<int, int, int, int>(minRow, maxRow, elutionStartCol, elutionEndCol);
        }

        private double GetSummedEnvelopeAtCharge(int row, int minCol, int maxCol, out int newMinCol, out int newMaxCol)
        {
            const int maxScanSkips = 2;
            const double goodEnoughBcDistance = 0.04;
            const double goodEnoughCorrCoeff = 0.7;
            
            var seedCol = -1;
            var seedDist = 10.0d;
            newMinCol = minCol;
            newMaxCol = maxCol;

            for (var j = minCol; j <= maxCol; j++)
            {
                if (!(_featureMatrix[row][j].AccurateMass > 0)) continue;
                var bcDist = _featureMatrix[row][j].DivergenceDist;

                if (bcDist < seedDist)
                {
                    seedCol = j;
                    seedDist = bcDist;
                }
            }

            if (seedCol < 0) return seedDist;

            // going forward
            var summedBcDist = 10.0d;
            var summedCorr = 0d;
            var tempEnvelope = new double[_theoreticalEnvelope.Size];
            var summedEnvelope = new double[_theoreticalEnvelope.Size];
            var n = 0;
            for(var col = seedCol; col < NColumns; col++)
            {
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
                    if (tempBcDist < summedBcDist)
                    {
                        summedBcDist = tempBcDist;
                        Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                    }
                    else
                    {
                        Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                        n++;
                    }                    
                }
            }
            return summedBcDist;
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
        }
    }
}
