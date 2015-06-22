using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security.Principal;
using System.Text;
using System.Threading.Tasks;
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
            

            _distProfileAcrossCharge = new double[NRows];
            _summedEnvelopeColRange = new int[NRows][];
            for(var i = 0; i < NRows; i++) _summedEnvelopeColRange[i] = new int[2];
            
            _xic = new double[NColumns];

            Comparer = new MzComparerWithBinning(NumOfBits);

            //_featureMatrixCreated = false;
            _featureMatrix = null;
            _refinedFeatureMatrix = null;
        }

        private const int NumOfBits = 27;

        public IList<LcMsPeakCluster> FindFeatures(int binNumber)
        {
            //var targetMass = Comparer.GetMzAverage(binNumber);
            //SetTargetMass(targetMass);
            var features = GetLcMs1PeakClusters(binNumber);

            return features;

            /*
            foreach (var feature in features)
            {
                FindLcMsPeakCluster(feature , (int)scan, (int)charge);    
            }*/

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

        public readonly int NColumns;
        public readonly int NRows;
        public readonly MzComparerWithBinning Comparer;

        private readonly int _maxThreadCount;
        private static List<Ms1Peak> _ms1PeakList;
        protected const int MaxEnvelopeLength = 30;

        private readonly double[] _distProfileAcrossCharge;
        private readonly int[][] _summedEnvelopeColRange;
        private readonly double[] _xic;

        //private bool _featureMatrixCreated;
        private LcMsPeakMatrixCell[][] _featureMatrix;
        private int[] _rows;
        private int[] _cols;

        private LcMsPeakMatrixCell[][] _refinedFeatureMatrix;
        //private double _refinedTargetMass;
        private int[] _refinedRows;
        private int[] _refinedCols;
        
        protected TheoreticalIsotopeEnvelope TheoreticalEnvelope;
        
        /*
        private bool _isAccurateMass;
        private void SetTargetMass(bool isAccurateMass, double targetMass, out int[] rows)
        {
            _isAccurateMass = isAccurateMass;
            var chargeLowerBound = Math.Floor(targetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Ceiling(targetMass / Run.MinMs1Mz);
            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);
            rows = Enumerable.Range((int)rowLb, (int)(rowUb - rowLb + 1)).ToArray();                

            if (_isAccurateMass)
            {
                _refinedTargetMass = targetMass;
                if (TheoreticalEnvelope == null) TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, 0.1d);
            }
            else
            {
                _targetMass = targetMass;
                TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, 0.1d);
            }
        }
        */
        private double _targetAccurateMass;

        private void BuildFeatureMatrix(bool isAccurateMass, double targetMass)
        {
            if (isAccurateMass)
            {
                if (_refinedFeatureMatrix == null)
                {
                    _refinedFeatureMatrix = new LcMsPeakMatrixCell[NRows][];
                    for (var i = 0; i < NRows; i++)
                    {
                        _refinedFeatureMatrix[i] = new LcMsPeakMatrixCell[NColumns];

                        for (var j = 0; j < NColumns; j++)
                            _refinedFeatureMatrix[i][j] = new LcMsPeakMatrixCell(MaxEnvelopeLength);
                    }
                }
                else
                {
                    // for quick computation
                    var smTol = new Tolerance(2.5);
                    var msTh = smTol.GetToleranceAsTh(targetMass);
                    if (Math.Abs(_targetAccurateMass - targetMass) < msTh) return;
                }
                _targetAccurateMass = targetMass;
            }
            else
            {
                if (_featureMatrix == null)
                {
                    _featureMatrix = new LcMsPeakMatrixCell[NRows][];
                    for (var i = 0; i < NRows; i++)
                    {
                        _featureMatrix[i] = new LcMsPeakMatrixCell[NColumns];
                        for (var j = 0; j < NColumns; j++)
                            _featureMatrix[i][j] = new LcMsPeakMatrixCell(MaxEnvelopeLength);
                    }
                }
            }

            //var tolerance = (isAccurateMass) ? new Tolerance(10) : new Tolerance(5);
            var tolerance = new Tolerance(5);
            var featureMatrix = (isAccurateMass) ? _refinedFeatureMatrix : _featureMatrix;

            var chargeLowerBound = Math.Floor(targetMass / Run.MaxMs1Mz);
            var chargeUpperBound = Math.Ceiling(targetMass / Run.MinMs1Mz);
            var rowLb = Math.Max(0, chargeLowerBound - MinScanCharge);
            var rowUb = Math.Min(NRows - 1, chargeUpperBound - MinScanCharge);
            var rows = Enumerable.Range((int)rowLb, (int)(rowUb - rowLb + 1)).ToArray();
            TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(targetMass, MaxEnvelopeLength, 0.1d);
            
            //SetTargetMass(isAccurateMass, targetMass, out rows);
            var queryMassBinNum = Comparer.GetBinNumber(targetMass);

            var options             = new ParallelOptions();
            var observedRows        = new bool[NRows];
            var observedCols        = new bool[NColumns];

            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            var mostAbuIsotopeInternalIndex = TheoreticalEnvelope.IndexOrderByRanking[0];
            var mostAbuIsotopeIndex = TheoreticalEnvelope.GetMostAbundantIsotope().Index;
            
            var minMs1Mz = _ms1PeakList.First().Mz;
            var maxMs1Mz = _ms1PeakList.Last().Mz;

            var signalToNoiseRatioCutoff = (isAccurateMass) ? 0 : 1.4826;
            var nPeaksCutoff = (targetMass > 2000) ? 3 : 2;

            Parallel.ForEach(rows, options, row =>
            {
                for (var col = 0; col < NColumns; col++) featureMatrix[row][col].Init();

                var charge = row + MinScanCharge;
                for (var k = 0; k < TheoreticalEnvelope.Size; k++)
                {
                    var i = TheoreticalEnvelope.IndexOrderByRanking[k];
                    var isotopeIndex = TheoreticalEnvelope.Isotopes[i].Index;

                    double isotopeMzLb;
                    double isotopeMzUb;

                    if (isAccurateMass)
                    {
                        var isotopeMz = Ion.GetIsotopeMz(targetMass, charge, isotopeIndex);
                        var mzTol = tolerance.GetToleranceAsTh(isotopeMz);
                        isotopeMzLb = isotopeMz - mzTol;
                        isotopeMzUb = isotopeMz + mzTol;
                    }
                    else
                    {
                        isotopeMzLb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum - 1), charge, isotopeIndex);
                        isotopeMzUb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum + 1), charge, isotopeIndex);
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
                            if (featureMatrix[row][col].EnvelopePeaks[i] == null || ms1Peak.Intensity > featureMatrix[row][col].EnvelopePeaks[i].Intensity)
                            {
                                featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                featureMatrix[row][col].AccurateMass = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                            }
                        }
                        else
                        {
                            if (isAccurateMass)
                            {
                                if (!(featureMatrix[row][col].AccurateMass > 0)) featureMatrix[row][col].AccurateMass = targetMass;
                            }
                            else
                            {
                                // mass binning search mode, most abundant peak should be exist
                                var mostAbuPeak = featureMatrix[row][col].EnvelopePeaks[mostAbuIsotopeInternalIndex];
                                if (mostAbuPeak == null) continue;
                                var expectedPeakMz = mostAbuPeak.Mz +
                                                     (Constants.C13MinusC12 * (isotopeIndex - mostAbuIsotopeIndex)) / charge;
                                if (Math.Abs(expectedPeakMz - ms1Peak.Mz) > tolerance.GetToleranceAsTh(ms1Peak.Mz))
                                    continue;                                
                            }

                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (featureMatrix[row][col].EnvelopePeaks[i] != null)
                            {
                                var tmpPeak = featureMatrix[row][col].EnvelopePeaks[i];
                                var bc1 = TheoreticalEnvelope.GetBhattacharyyaDistance(featureMatrix[row][col].EnvelopePeaks);
                                featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;
                                var bc2 = TheoreticalEnvelope.GetBhattacharyyaDistance(featureMatrix[row][col].EnvelopePeaks);
                                if (bc1 < bc2) featureMatrix[row][col].EnvelopePeaks[i] = tmpPeak;
                            }
                            else
                            {
                                featureMatrix[row][col].EnvelopePeaks[i] = ms1Peak;

                            }
                        }
                    }
                }

                for (var col = 0; col < NColumns; col++)
                {
                    if (!(featureMatrix[row][col].AccurateMass > 0)) continue;

                    var nPeaks = featureMatrix[row][col].EnvelopePeaks.Count(p => p != null && p.Active);
                    if (isAccurateMass)
                    {
                        if (nPeaks >= nPeaksCutoff)
                        {
                            featureMatrix[row][col].CorrelationCoeff = TheoreticalEnvelope.GetPearsonCorrelation(featureMatrix[row][col].EnvelopePeaks);
                            featureMatrix[row][col].DivergenceDist = TheoreticalEnvelope.GetBhattacharyyaDistance(featureMatrix[row][col].EnvelopePeaks);
                            if (!observedRows[row]) observedRows[row] = true;
                            if (!observedCols[col]) observedCols[col] = true;
                        }
                        else
                        {
                            featureMatrix[row][col].AccurateMass = 0d;
                        }
                    }
                    else
                    {
                        var mostAbuPeakIntensity = featureMatrix[row][col].EnvelopePeaks[mostAbuIsotopeInternalIndex].Intensity;
                        var signalToNoiseRatio = mostAbuPeakIntensity / Ms1Spectra[col].MedianIntensity;
                        if (signalToNoiseRatio > signalToNoiseRatioCutoff && nPeaks >= nPeaksCutoff)
                        {
                            featureMatrix[row][col].CorrelationCoeff = TheoreticalEnvelope.GetPearsonCorrelation(featureMatrix[row][col].EnvelopePeaks);
                            featureMatrix[row][col].DivergenceDist = TheoreticalEnvelope.GetBhattacharyyaDistance(featureMatrix[row][col].EnvelopePeaks);
                            if (!observedRows[row]) observedRows[row] = true;
                            if (!observedCols[col]) observedCols[col] = true;
                        }
                        else
                        {
                            featureMatrix[row][col].AccurateMass = 0d;
                        }
                    }
                }
            }// end or row for-loop
            );

            var temp = new List<int>();
            for (var i = 0; i < observedRows.Length; i++) if (observedRows[i]) temp.Add(i);
            
            if (isAccurateMass) _refinedRows = temp.ToArray();
            else _rows = temp.ToArray();

            temp.Clear();
            for (var i = 0; i < observedCols.Length; i++) if (observedCols[i]) temp.Add(i);

            if (isAccurateMass) _refinedCols = temp.ToArray();
            else _cols = temp.ToArray();
        }

        private IEnumerable<ObservedIsotopeEnvelope> GetSeedCells(double targetMass)
        {
            var seedList = new List<KeyValuePair<double, ObservedIsotopeEnvelope>>();
            var ms1ScanNums = Run.GetMs1ScanVector();
            //var mostAbuIsotopeInternalIndex = TheoreticalEnvelope.IndexOrderByRanking[0];

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

                    //var mostAbuPeakIntensity = _featureMatrix[i][j].EnvelopePeaks[mostAbuIsotopeInternalIndex].Intensity;
                    //var signalToNoiseRatio = mostAbuPeakIntensity / Ms1Spectra[j].MedianIntensity;
                    //if (signalToNoiseRatio < 3) continue;

                    if (_featureMatrix[i][j].EnvelopePeaks.Count(p => p != null && p.Active) < nPeaksCutoff) continue;

                    var bcDist = _featureMatrix[i][j].DivergenceDist;
                    var corr = _featureMatrix[i][j].CorrelationCoeff;

                    if (bcDist > 0.3 && corr < 0.4) continue;

                    var seed = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, TheoreticalEnvelope);
                    seedList.Add(new KeyValuePair<double, ObservedIsotopeEnvelope>(bcDist, seed));
                }
            }
            return seedList.OrderBy(x => x.Key).Select(x => x.Value);
        }

        private IList<LcMsPeakCluster> GetLcMs1PeakClusters(int binNumber)
        {
            var targetMass = Comparer.GetMzAverage(binNumber);
            BuildFeatureMatrix(false, targetMass); // should be called first

            var clusters = new List<LcMsPeakCluster>();
            var tempEnvelope = new double[TheoreticalEnvelope.Size];
            var tempEnvelope2 = new double[TheoreticalEnvelope.Size];
            
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
                    var minRw = Math.Max(charge - MinScanCharge - chargeNeighborGap, _rows.First());
                    var maxRw = Math.Min(charge - MinScanCharge + chargeNeighborGap, _rows.Last());

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

                            //if (_featureMatrix[i][j].DivergenceDist > GetBcDistTh(newCluster.NetLength) || _featureMatrix[i][j].CorrelationCoeff < GetCorrTh(newCluster.NetLength)) continue;
                            // Summing envelope from paired isotopic envelopes 
                            
                            //cell.Peaks.SumEnvelopeTo(tempEnvelope);
                            //_featureMatrix[i][j].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);

                            Array.Copy(tempEnvelope, tempEnvelope2, tempEnvelope2.Length);

                            _featureMatrix[i][j].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                            var newDivergence = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                            var newCorrelation = TheoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);


                            if (_featureMatrix[i][j].DivergenceDist < 0.04 &&
                                _featureMatrix[i][j].CorrelationCoeff > 0.9 || newDivergence < summedBcDist ||
                                newCorrelation > summedCorr)
                            {
                                var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass,
                                    i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks,
                                    TheoreticalEnvelope);

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
                    var refinedCluster = GetLcMsPeakCluster(newCluster.RepresentativeMass,
                        newCluster.RepresentativeCharge, newCluster.MinScanNum, newCluster.MaxScanNum);
                    if (refinedCluster != null)
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

                //for (var i = newCluster.MinCharge - MinScanCharge; i <= newCluster.MaxCharge - MinScanCharge; i++)
                    //for (var j = ms1ScanNumToIndex[newCluster.MinScanNum]; j <= ms1ScanNumToIndex[newCluster.MaxScanNum]; j++) _featureMatrix[i][j].CheckedOutFlag = true;
            }
            return clusters;            
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

            for (var j = col - 1; j >= _refinedCols.First(); j--)
            {
                if (elutionTime - Run.GetElutionTime(ms1ScanNums[j]) > halfWindowSize && (col - j) >= minScanCount) break;
                if (j < minCol) minCol = j;
            }

            for (var j = col + 1; j < _refinedCols.Last(); j++)
            {
                if (Run.GetElutionTime(ms1ScanNums[j]) - elutionTime > halfWindowSize && (j - col) >= minScanCount) break;
                if (j > maxCol) maxCol = j;
            }
            return new Tuple<int, int>(minCol, maxCol);
        }
        
        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);
        
        //private const double InitialElutionPeriodUb = 5.0;
        //private const double InitialElutionPeriodLb = 0.5;
        public LcMsPeakCluster GetLcMsPeakCluster(double targetMass, int targetCharge, int targetMinScanNum, int targetMaxScanNum, int seedScanNum = -1)
        {
            BuildFeatureMatrix(true, targetMass); // should be called first

            //set initial elution period for initial sampling
          
            var ms1ScanNums = Run.GetMs1ScanVector();
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            
            if (Run.GetMsLevel(targetMinScanNum) > 1) targetMinScanNum = Run.GetPrevScanNum(targetMinScanNum, 1);
            if (Run.GetMsLevel(targetMaxScanNum) > 1) targetMaxScanNum = Run.GetNextScanNum(targetMaxScanNum, 1);

            var seedCol = -1;
            var row = targetCharge - MinScanCharge;

            if (seedScanNum < 0)
            {
                var targetScanNum = (int)((targetMinScanNum + targetMaxScanNum) * 0.5d);
                
                seedCol = Array.BinarySearch(ms1ScanNums, targetScanNum);
                if (seedCol < 0) seedCol = ~seedCol;
            }
            else
            {
                seedCol = ms1ScanNumToIndex[seedScanNum];
            }

            var initialElutionPeriod = Math.Max(Math.Min(Run.GetElutionTime(Run.MaxLcScan) * 0.003, 5.0), 0.5);
            var tempWindow = GetElutionWindow(seedCol, initialElutionPeriod);
            var minCol = Math.Min(tempWindow.Item1, ms1ScanNumToIndex[targetMinScanNum]);
            var maxCol = Math.Max(tempWindow.Item2, ms1ScanNumToIndex[targetMaxScanNum]);
            
            // 1) determine charge state range and elution period
            var range = DetermineFeatureRange(row, minCol, maxCol, seedCol);
            var minRow = range.Item1;
            var maxRow = range.Item2;
            minCol = range.Item3;
            maxCol = range.Item4;
            
            // 3) collect envelopes
            var feature = CollectLcMsPeaks(minRow, maxRow, minCol, maxCol);

            if (feature == null) return null;

            feature.UpdateScore(Ms1Spectra);

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

        public double[] GetEntireXic(LcMsPeakCluster feature)
        {
            Array.Clear(_xic, 0, _xic.Length);
            var row = feature.RepresentativeCharge - MinScanCharge;
            foreach (var col in _refinedCols)
            {
                _xic[col] = _refinedFeatureMatrix[row][col].Intensity;
            }
            var smoothedXic = Smoother.Smooth(_xic);
            return smoothedXic;
        }

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

            feature.EnvelopeDistanceScoreAcrossCharge = new double[maxRow - minRow + 1];
            feature.EnvelopeIntensityScoreAcrossCharge = new double[maxRow - minRow + 1];
            feature.NormalizedIntensityAcrossCharge = new double[maxRow - minRow + 1];
            
            var intensities = new double[TheoreticalEnvelope.Size];

            var comparer = new MzComparerWithBinning(28);
            var NumMzBins = comparer.GetBinNumber(Run.MaxMs1Mz) - comparer.GetBinNumber(Run.MinMs1Mz) + 1;
            
            for(var i = minRow; i <= maxRow; i++)
            {
                Array.Clear(intensities, 0, intensities.Length);
                Array.Clear(summedPeakCount, 0, summedPeakCount.Length);
                Array.Clear(summedPeakIntensity, 0, summedPeakIntensity.Length);
                Array.Clear(summedMedianIntensity, 0, summedMedianIntensity.Length);

                for (var j = minCol; j <= maxCol; j++)
                {
                    for (var k = 0; k < TheoreticalEnvelope.Size; k++)
                    {
                        var r = rnd.Next(0, NumMzBins);
                        if (r < Ms1Spectra[j].Peaks.Length)
                        {
                            r = rnd.Next(0, n);
                            var randomPeak = _ms1PeakList[r];
                            intensities[k] += randomPeak.Intensity;
                            
                            if (k == TheoreticalEnvelope.IndexOrderByRanking[0])
                            {
                                summedPeakCount[i]++;
                                summedPeakIntensity[i] += randomPeak.Intensity;
                                summedMedianIntensity[i] += Ms1Spectra[randomPeak.Ms1SpecIndex].GetLocalMedianIntensity(randomPeak, feature.Mass);
                            }
                        }
                        
                    }
                }
                
                feature.EnvelopeDistanceScoreAcrossCharge[i - minRow] = TheoreticalEnvelope.GetBhattacharyyaDistance(intensities);
                if (summedMedianIntensity[i] > 0) feature.EnvelopeIntensityScoreAcrossCharge[i - minRow] = summedPeakIntensity[i] / summedMedianIntensity[i];
                feature.NormalizedIntensityAcrossCharge[i - minRow] = intensities.Sum();
            }

            var s = feature.NormalizedIntensityAcrossCharge.Sum();
            foreach (var i in _refinedRows)
                feature.NormalizedIntensityAcrossCharge[i - minRow] = feature.NormalizedIntensityAcrossCharge[i - minRow] / s;

            feature.SetChargeRange(minRow + MinScanCharge, maxRow + MinScanCharge); 
        }
      

        private LcMsPeakCluster CollectLcMsPeaks(int minRow, int maxRow, int minCol, int maxCol)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();

            var envelopes = new List<ObservedIsotopeEnvelope>();
            var bestBcDist = 100d;
            ObservedIsotopeEnvelope bestEnvelope = null;

            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!_refinedFeatureMatrix[i][j].Exist) continue;

                    var envelope = new ObservedIsotopeEnvelope(_refinedFeatureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _refinedFeatureMatrix[i][j].EnvelopePeaks, TheoreticalEnvelope);
                    envelopes.Add(envelope);

                    if (_refinedFeatureMatrix[i][j].DivergenceDist < bestBcDist)
                    {
                        bestBcDist = _refinedFeatureMatrix[i][j].DivergenceDist;
                        bestEnvelope = envelope;
                    }
                }
            }

            if (bestEnvelope == null)
            {
                return null;
            }

            var cluster = new LcMsPeakCluster(Run, bestEnvelope);
            foreach(var e in envelopes) cluster.AddObservedEnvelope(e);
            
            return cluster;
        }

        private Tuple<int, int, int, int> DetermineFeatureRange(int seedRow, int seedMinCol, int seedMaxCol, int seedCol)
        {
            //const double searchStopBcThreshold = 0.3;
            var bestBcDist = 10.0d;
            var bestBcDistRow = 0;
            Array.Clear(_distProfileAcrossCharge, 0, _distProfileAcrossCharge.Length);
            foreach (var i in _refinedRows)
            {
                var c1 = 0;
                var c2 = 0;
                var summedBcDist = GetSummedEnvelopeAtCharge(i, seedMinCol, seedMaxCol, out c1, out c2);
                _distProfileAcrossCharge[i] = summedBcDist;

                if (c1 <= seedCol && seedCol <= c2 && _distProfileAcrossCharge[i] < bestBcDist)
                {
                    bestBcDistRow = i;
                    bestBcDist = _distProfileAcrossCharge[i];
                }

                _summedEnvelopeColRange[i][0] = c1;
                _summedEnvelopeColRange[i][1] = c2;
            }

            //var bcCutoff = OtsuThreshold.GetThreshold(_distProfileAcrossCharge, 0, searchStopBcThreshold, 0.0025, 1e-9, searchStopBcThreshold);
            //bcCutoff = Math.Min(Math.Max(bcCutoff, 0.05), 0.15);
            const double bcCutoff = 0.1;
            var minCol = NColumns - 1;
            var maxCol = 0;

            var minColBestCharge = _summedEnvelopeColRange[bestBcDistRow][0];
            var maxColBestCharge = _summedEnvelopeColRange[bestBcDistRow][1];

            // determine temporary charge and scan boundaries
            var minRow = seedRow;
            var maxRow = seedRow;
            foreach (var i in _refinedRows)
            {
                if (_distProfileAcrossCharge[i] > bcCutoff) continue;

                if (_summedEnvelopeColRange[i][1] < minColBestCharge) continue;
                if (_summedEnvelopeColRange[i][0] > maxColBestCharge) continue;

                if (i < minRow) minRow = i;
                if (i > maxRow) maxRow = i;
                
                if (_summedEnvelopeColRange[i][0] < minCol) minCol = _summedEnvelopeColRange[i][0];
                if (_summedEnvelopeColRange[i][1] > maxCol) maxCol = _summedEnvelopeColRange[i][1];
            }

            // determine Min and Max Columns (based on elution profile)
            Array.Clear(_xic, 0, _xic.Length);
            for (var j = minCol; j <= maxCol; j++)
            {
                _xic[j] = _refinedFeatureMatrix[bestBcDistRow][j].Intensity;
            }
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

            var ms1ScanNums = Run.GetMs1ScanVector();
            var oneSigIntensity = apexIntensity * Math.Exp(-1);
            //var twoSigIntensity = apexIntensity * Math.Exp(-2);
            var threeSigIntensity = apexIntensity * Math.Exp(-4.5);


            // going backward in elution time
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

            // find two sigma point at elution time axis
            var elutionStartColByTwoSigma = elutionStartColByOneSigma;
            for (var j = elutionStartColByOneSigma - 1; j >= 0; j--)
            {
                if (Run.GetElutionTime(ms1ScanNums[j]) < Run.GetElutionTime(ms1ScanNums[apexCol]) - 2*oneSigPeriod)
                    elutionStartColByTwoSigma = j;
                else break;
            }

            // possible extension?
            var elutionStartCol = elutionStartColByTwoSigma;
            for (var j = elutionStartColByTwoSigma - 1; j >= minCol; j--)
            {
                if (smoothedXic[j] < threeSigIntensity) break;

                var needBreak = true;
                for (var k = j; k >= j - 3; k--)
                {
                    if (k < minCol) break;
                    if (smoothedXic[k] < smoothedXic[j + 1]) needBreak = false;
                }
                if (needBreak) break;
                    
                elutionStartCol = j;
            }

            // going forward in elution time
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
            const double goodEnoughBcDistance = 0.02;
            const double goodEnoughCorrCoeff = 0.7;
            
            var seedCol = -1;
            var seedDist = 10.0d;
            newMinCol = minCol;
            newMaxCol = maxCol;

            for (var j = minCol; j <= maxCol; j++)
            {
                if (!(_refinedFeatureMatrix[row][j].AccurateMass > 0)) continue;
                var bcDist = _refinedFeatureMatrix[row][j].DivergenceDist;

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
            var tempEnvelope = new double[TheoreticalEnvelope.Size];
            var summedEnvelope = new double[TheoreticalEnvelope.Size];
            var n = 0;
            for(var col = seedCol; col < NColumns; col++)
            {
                if (n >= maxScanSkips)
                {
                    newMaxCol = col;
                    break;
                }

                if (!_refinedFeatureMatrix[row][col].Exist)
                {
                    n++;
                    continue;
                }

                _refinedFeatureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                var tempBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                var tempCorr = TheoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);

                if (tempBcDist < summedBcDist || tempCorr > summedCorr || _refinedFeatureMatrix[row][col].DivergenceDist < goodEnoughBcDistance || _refinedFeatureMatrix[row][col].CorrelationCoeff > goodEnoughCorrCoeff)
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

                if (!_refinedFeatureMatrix[row][col].Exist)
                {
                    n++;
                    continue;
                }

                _refinedFeatureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                var tempBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                var tempCorr = TheoreticalEnvelope.GetPearsonCorrelation(tempEnvelope);

                if (tempBcDist < summedBcDist || tempCorr > summedCorr || _refinedFeatureMatrix[row][col].DivergenceDist < goodEnoughBcDistance || _refinedFeatureMatrix[row][col].CorrelationCoeff > goodEnoughCorrCoeff)
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
            _refinedFeatureMatrix[row][seedCol].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
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

                    if (!_refinedFeatureMatrix[row][col].Exist) continue;

                    _refinedFeatureMatrix[row][col].EnvelopePeaks.SumEnvelopeTo(tempEnvelope);
                    var tempBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
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
        /*
      private void EvaludateSummedEnvelope(LcMsPeakCluster feature)
      {
          var ms1ScanToIndex = Run.GetMs1ScanNumToIndex();
          var ms1ScanVector = Run.GetMs1ScanVector();
            
          var minRow = feature.MinCharge - MinScanCharge;
          var maxRow = feature.MaxCharge - MinScanCharge;
            
          var minCol = ms1ScanToIndex[feature.MinScanNum];
          var maxCol = ms1ScanToIndex[feature.MaxScanNum];

          var summedPeakCount         = new int[NRows];
          var summedPeakIntensity     = new double[NRows];
          var summedMedianIntensity   = new double[NRows];

          var mostAbuIdx = feature.Envelopes[0].TheoreticalEnvelope.IndexOrderByRanking[0];

          feature.EnvelopeDistanceScoreAcrossCharge = new double[maxRow - minRow + 1];
          feature.EnvelopeIntensityScoreAcrossCharge = new double[maxRow - minRow + 1];
          feature.NormalizedIntensityAcrossCharge = new double[maxRow - minRow + 1];
            
          var intensities = new double[TheoreticalEnvelope.Size];

          var tempBestBcDistEnvelope = new double[TheoreticalEnvelope.Size];
          var tempBestBcDist = 10.0d;
          var tempBestBcDistRow = 0;

          for (var i = minRow; i <= maxRow; i++)
          {
              Array.Clear(intensities, 0, intensities.Length);
              Array.Clear(summedPeakCount, 0, summedPeakCount.Length);
              Array.Clear(summedPeakIntensity, 0, summedPeakIntensity.Length);
              Array.Clear(summedMedianIntensity, 0, summedMedianIntensity.Length);

              var bestBcDist = 10.0d;
              for (var j = minCol; j <= maxCol; j++)
              {
                  if (!_featureMatrix[i][j].Exist) continue;
                    
                  _featureMatrix[i][j].EnvelopePeaks.SumEnvelopeTo(intensities);
                  var mostAbuPeak = _featureMatrix[i][j].EnvelopePeaks[mostAbuIdx];

                  var newBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(_featureMatrix[i][j].EnvelopePeaks);
                  if (newBcDist < bestBcDist) bestBcDist = newBcDist;

                  if (mostAbuPeak != null && mostAbuPeak.Active)
                  {
                      summedPeakCount[i]++;
                      summedPeakIntensity[i] += mostAbuPeak.Intensity;
                      summedMedianIntensity[i] +=
                          Ms1Spectra[mostAbuPeak.Ms1SpecIndex].GetLocalMedianIntensity(mostAbuPeak, TargetMass);                        
                  }
              }

              feature.EnvelopeDistanceScoreAcrossCharge[i - minRow] = Math.Min(TheoreticalEnvelope.GetBhattacharyyaDistance(intensities), bestBcDist);
              if (feature.EnvelopeDistanceScoreAcrossCharge[i - minRow] < tempBestBcDist)
              {
                  tempBestBcDist = feature.EnvelopeDistanceScoreAcrossCharge[i - minRow];
                  tempBestBcDistRow = i;
                  Array.Copy(intensities, tempBestBcDistEnvelope, tempBestBcDistEnvelope.Length);
              }

              if (summedMedianIntensity[i] > 0)
                  feature.EnvelopeIntensityScoreAcrossCharge[i - minRow] = summedPeakIntensity[i] / summedMedianIntensity[i];
              feature.NormalizedIntensityAcrossCharge[i - minRow] = intensities.Sum();
          }
            
          var s = feature.NormalizedIntensityAcrossCharge.Sum();
          for (var i = minRow; i <= maxRow; i++)
          {
              feature.NormalizedIntensityAcrossCharge[i - minRow] = feature.NormalizedIntensityAcrossCharge[i - minRow] / s;
          }

          tempBestBcDist = 10.0d;
          var tempBestBcDistCol = 0;
          for (var j = minCol; j < maxCol; j++)
          {
              if (_featureMatrix[tempBestBcDistRow][j].DivergenceDist < tempBestBcDist)
              {
                  tempBestBcDist = _featureMatrix[tempBestBcDistRow][j].DivergenceDist;
                  tempBestBcDistCol = j;
              }
          }

          var tempBestBcDistMz = 0d;
          var peaks = _featureMatrix[tempBestBcDistRow][tempBestBcDistCol].EnvelopePeaks;
          foreach (var i in TheoreticalEnvelope.IndexOrderByRanking)
          {
              if (peaks[i] != null && peaks[i].Active)
              {
                  tempBestBcDistMz = peaks[i].Mz;
              }
          }
            
          feature.SetRepresentativeEnvelope(tempBestBcDistRow + MinScanCharge, tempBestBcDistMz, ms1ScanVector[tempBestBcDistCol], tempBestBcDistEnvelope);
            
      }*/
        /*
               private LcMsPeakCluster CollectLcMsPeaks_Selective(int minRow, int maxRow, int minCol, int maxCol)
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
                           if (!(_featureMatrix[i][j].AccurateMass > 0)) continue;
                           var envelope = new ObservedIsotopeEnvelope(_featureMatrix[i][j].AccurateMass, i + MinScanCharge, ms1ScanNums[j], _featureMatrix[i][j].EnvelopePeaks, TheoreticalEnvelope);
                           envelopes.Add(envelope);
                           bcDistList.Add(_featureMatrix[i][j].DivergenceDist);
                           corrList.Add(_featureMatrix[i][j].CorrelationCoeff);

                           if (_featureMatrix[i][j].DivergenceDist < bestBcDist)
                           {
                               bestBcDist = _featureMatrix[i][j].DivergenceDist;
                               bestEnvelope = envelope;
                           }
                       }
                   }

                   if (bestEnvelope == null)
                   {
                       return null;
                   }
            
                   var bcDistSample = bcDistList.Where(d => d < bcDistList.Median()).ToList();
                   var meanStd = bcDistSample.MeanStandardDeviation();
                   var bcCutoff = Math.Max(meanStd.Item1 + 2.5 * meanStd.Item2, 0.04);
                   bcCutoff = Math.Min(bcCutoff, 0.5);

                   var corrSample = corrList.Where(d => d < corrList.Median()).ToList();
                   meanStd = corrSample.MeanStandardDeviation();
                   var corrCutoff = Math.Max(meanStd.Item1 + 2.5 * meanStd.Item2, 0.2);
                   corrCutoff = Math.Min(bcCutoff, 0.9);

                   var cluster = new LcMsPeakCluster(Run, bestEnvelope);
                   for (var k = 0; k < bcDistList.Count; k++)
                   {
                       if (bcDistList[k] < bcCutoff || corrList[k] > corrCutoff) cluster.AddObservedEnvelope(envelopes[k]); 
                   }

                   return cluster;
               }
               */
        
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

            internal bool Exist { get { return (AccurateMass > 0d); } }

            internal double Intensity
            {
                get
                {
                    return EnvelopePeaks.Where(envelopePeak => envelopePeak != null && envelopePeak.Active).Sum(envelopePeak => envelopePeak.Intensity);
                }
            }
        }
    }
}
