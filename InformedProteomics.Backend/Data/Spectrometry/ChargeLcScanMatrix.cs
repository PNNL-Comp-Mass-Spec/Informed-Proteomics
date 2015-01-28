using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter, IMs1FeatureExtract
    {
        public ChargeLcScanMatrix(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = -1, int numBits = 27)
        {
            _run = run;
            _minScanCharge = minScanCharge;
            _maxScanCharge = maxScanCharge;
            _maxThreadCount = maxThreadCount;
            //_predictor      = null;

            _comparer = new MzComparerWithBinning(numBits);
            _ms1ScanNums = run.GetMs1ScanVector();
            _nScans = _ms1ScanNums.Length;

            var ms1PeakList = new List<Peak>();
            var ms1PeakScanIndex = new List<Ms1PeakInfo>();
            _cachedSpectrum = new List<ChargeLcScanSpectrum>();
            for (var i = 0; i < _ms1ScanNums.Length; i++)
            {
                var spec = new ChargeLcScanSpectrum(_run.GetSpectrum(_ms1ScanNums[i]));
                _cachedSpectrum.Add(spec);
                ms1PeakList.AddRange(spec.Peaks);
                for (var j = 0; j < spec.Peaks.Length; j++) ms1PeakScanIndex.Add(new Ms1PeakInfo(i, j));
            }

            _ms1PeakArr = ms1PeakList.ToArray();
            _ms1PeakInfoArr = ms1PeakScanIndex.ToArray();
            Array.Sort(_ms1PeakArr, _ms1PeakInfoArr);

            _mostAbuIsotopePeakIndex = new int[MaxChargeLength][];
            _correlationMap = new double[MaxChargeLength][];
            _featureMatrix = new double[MaxChargeLength][][];
            _checkedOut = new bool[MaxChargeLength][];
            _accurateMass = new double[MaxChargeLength][];

            _highestPeakIndex = new int[MaxChargeLength][];
            _highestPeakIntensity = new double[MaxChargeLength][];

            for (var i = 0; i < MaxChargeLength; i++)
            {
                _mostAbuIsotopePeakIndex[i] = new int[_nScans];
                _checkedOut[i] = new bool[_nScans];
                _correlationMap[i] = new double[_nScans];
                _featureMatrix[i] = new double[_nScans][];
                _highestPeakIndex[i] = new int[_nScans];
                _highestPeakIntensity[i] = new double[_nScans];
                _accurateMass[i] = new double[_nScans];

                for (var j = 0; j < _nScans; j++)
                    _featureMatrix[i][j] = new double[MaxEnvelopeLength];
            }
        }

        public string GetFeatureFile(string rawFilePath, double minSearchMass = 3000, double maxSearchMass = 50000)
        {
            var outTsvFilePath = GetFeatureFilePath(rawFilePath);
            if (File.Exists(outTsvFilePath)) return outTsvFilePath;
            return GenerateFeatureFile(rawFilePath, minSearchMass, maxSearchMass, false, 0.2, false, true);
        }

        public string GenerateFeatureFile(string rawFilePath, double minMass = 3000, double maxMass = 50000, bool massCollapse = false, 
            double probabilityThreshold = 0.5, bool scoreReport = false, bool csvOutput = false, bool oldFormat = false)
        {
            var j = rawFilePath.LastIndexOf('.');
            var outPath = rawFilePath.Substring(0, j);

            var outTsvFilePath = GetFeatureFilePath(rawFilePath);
            var outCsvFilePath = string.Format("{0}_ms1ft.csv", outPath);
            
            _minSearchMassBin = _comparer.GetBinNumber(minMass);
            _maxSearchMassBin = _comparer.GetBinNumber(maxMass);
            //_predictor = predictor;

            _massBinToClusterMap = new List<ChargeLcScanCluster>[_maxSearchMassBin - _minSearchMassBin + 1];
            for (var i = 0; i < _massBinToClusterMap.Length; i++) _massBinToClusterMap[i] = new List<ChargeLcScanCluster>();

            StreamWriter csvWriter = null;
            if (csvOutput)
            {
                csvWriter = new System.IO.StreamWriter(outCsvFilePath);
                csvWriter.WriteLine("scan_num,charge,monoisotopic_mw,mz,abundance,probability");
            }

            var tsvWriter = new System.IO.StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(ChargeLcScanCluster.GetHeaderString(scoreReport, oldFormat));

            Console.WriteLine("Start MS1 feature extracting...");
            Console.WriteLine("Mass Range {0} ~ {1}", minMass, maxMass);
            Console.WriteLine("Charge Range {0} ~ {1}", _minScanCharge, _maxScanCharge);
            Console.WriteLine("Output File\t{0}", outTsvFilePath);

            if (csvOutput) Console.WriteLine("Csv Output File\t{0}", outCsvFilePath);

            var binNumCursor = _minSearchMassBin;
            var stopwatch = Stopwatch.StartNew();
            stopwatch.Stop(); stopwatch.Reset();

            stopwatch.Start();
            _nOutFatures = 0;
            for (var binNum = _minSearchMassBin; binNum <= _maxSearchMassBin; binNum++)
            {
                ExtractFeatures(binNum, massCollapse, probabilityThreshold);

                if (binNum > _minSearchMassBin && (binNum - _minSearchMassBin) % 3000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - _minSearchMassBin;
                    //var remaining = (totalMassBins - processedBins) * (elapsed / processedBins);
                    FlushOutput(binNumCursor, binNum - 500, tsvWriter, probabilityThreshold, scoreReport, csvWriter);
                    Console.WriteLine("Processing {0} mass bins({1:0.00} Da). Elapsed Time = {2:0.000} sec. # of features = {3}", processedBins, _comparer.GetMzEnd(binNum), elapsed, _nOutFatures);
                    binNumCursor = binNum - 500 + 1;
                }
            }

            FlushOutput(binNumCursor, _maxSearchMassBin, tsvWriter, probabilityThreshold, scoreReport, csvWriter);
            tsvWriter.Close();
            if (csvOutput) csvWriter.Close();

            Console.WriteLine("Complete MS1 feature finding...Total Extracted Features = {0}...Elapsed Time = {1:0.000} [sec]", _nOutFatures, (stopwatch.ElapsedMilliseconds) / 1000.0d);

            return outTsvFilePath;
            //return _nOutFatures;
        }

        public IEnumerable<ChargeLcScanCluster> GetProbableChargeScanClusters(int queryMassBinNum)
        {
            return GetProbableChargeScanClusters(_comparer.GetMzAverage(queryMassBinNum));
        }

        public IEnumerable<Ms1Feature> GetProbableChargeScanRegions(double queryMass)
        {
            return GetProbableChargeScanClusters(queryMass);
        }

        public IEnumerable<ChargeLcScanCluster> GetProbableChargeScanClusters(double queryMass)
        {
            SetQueryMass(queryMass);
            var clusters = FindClusters();
            clusters = GetFilteredClusters(clusters);
            return clusters;
        }

        public IEnumerable<ChargeLcScanCluster> GetAllClusters(double queryMass)
        {
            SetQueryMass(queryMass);
            var clusters = FindClusters();
            clusters = GetFilteredClusters(clusters, false);
            return clusters;
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(ChargeLcScanCluster cluster)
        {
            var ms2ScanNumSet = new List<int>();
            for (var charge = _minScanCharge; charge <= _maxScanCharge; charge++)
            {
                var mostAbuMz = Ion.GetIsotopeMz(cluster.RepresentativeMass, charge, _isotopeList.GetMostAbundantIsotope().Index);
                var ms2ScanNums = _run.GetFragmentationSpectraScanNums(mostAbuMz);
                ms2ScanNumSet.AddRange(ms2ScanNums.Where(sn => sn > cluster.MinScanNum && sn < cluster.MaxScanNum));
            }

            return ms2ScanNumSet.Distinct();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double queryMass)
        {
            var clusters = GetProbableChargeScanClusters(queryMass);
            var ms2ScanNumSet = new List<int>();
            foreach (var cluster in clusters)
            {
                for (var charge = cluster.MinCharge; charge <= cluster.MaxCharge; charge++)
                {
                    var mostAbuMz = Ion.GetIsotopeMz(queryMass, charge, _isotopeList.GetMostAbundantIsotope().Index);
                    var ms2ScanNums = _run.GetFragmentationSpectraScanNums(mostAbuMz);
                    ms2ScanNumSet.AddRange(ms2ScanNums.Where(sn => sn > cluster.MinScanNum && sn < cluster.MaxScanNum));
                }
            }
            return ms2ScanNumSet.Distinct();
        }

        public static MzComparerWithBinning GetMzComparerWithBinning()
        {
            return _comparer;
        }

        private static LcMsRun _run;
        private static int _minScanCharge;
        private static int _maxScanCharge;
        private static int _nScans;
        private static int[] _ms1ScanNums;
        private static int _maxThreadCount;

        private static MzComparerWithBinning _comparer;

        internal static Tolerance MzTolerance = new Tolerance(10);
        private static List<ChargeLcScanSpectrum> _cachedSpectrum;

        private static Peak[] _ms1PeakArr;
        private static Ms1PeakInfo[] _ms1PeakInfoArr;

        private const int MaxEnvelopeLength = 30;
        private const int MaxChargeLength = 40;

        private readonly int[][] _mostAbuIsotopePeakIndex;
        private readonly double[][] _correlationMap;
        private readonly double[][][] _featureMatrix;
        private readonly bool[][] _checkedOut;

        private readonly double[][] _accurateMass;

        private readonly double[][] _highestPeakIntensity;
        private readonly int[][] _highestPeakIndex;

        private IsotopeList _isotopeList;
        private int[] _chargeIndexes;
        private IntRange _chargeRange;
        private double _queryMass;

        private string GetFeatureFilePath(string rawFilePath)
        {
            var j = rawFilePath.LastIndexOf('.');
            var outPath = rawFilePath.Substring(0, j);
            var outTsvFilePath = string.Format("{0}.ms1ft", outPath);
            return outTsvFilePath;
        }

        private void SetQueryMass(double queryMass)
        {
            _queryMass = queryMass;
            _chargeRange = GetScanChargeRange(_queryMass);
        }

        private IntRange GetScanChargeRange(double mass)
        {
            if (mass < 5000.0d) return new IntRange(_minScanCharge, 10);

            var chargeLb = (int)Math.Max(_minScanCharge, Math.Floor((13.0 / 2.5) * (mass / 10000d) - 0.6));
            var chargeUb = (int)Math.Min(_maxScanCharge, Math.Ceiling(18 * (mass / 10000d) + 8));

            if (chargeUb - chargeLb + 1 > MaxChargeLength) chargeUb = chargeLb + MaxChargeLength - 1;

            return new IntRange(chargeLb, chargeUb);
        }

        private void BuildFeatureMatrix()
        {
            var queryMassBinNum = _comparer.GetBinNumber(_queryMass);
            var observedCharges = new bool[_chargeRange.Length];

            _isotopeList = new IsotopeList(_comparer.GetMzAverage(queryMassBinNum), MaxEnvelopeLength);

            var rows = Enumerable.Range(0, _chargeRange.Length);

            var options = new ParallelOptions();
            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            //foreach(var row in rows)
            Parallel.ForEach(rows, options, row =>
            {
                Array.Clear(_mostAbuIsotopePeakIndex[row], 0, _nScans);
                Array.Clear(_correlationMap[row], 0, _nScans);
                Array.Clear(_checkedOut[row], 0, _nScans);

                Array.Clear(_highestPeakIndex[row], 0, _nScans);
                Array.Clear(_highestPeakIntensity[row], 0, _nScans);
                Array.Clear(_accurateMass[row], 0, _nScans);

                for (var col = 0; col < _nScans; col++)
                {
                    Array.Clear(_featureMatrix[row][col], 0, _featureMatrix[row][col].Length);
                }

                var charge = row + _chargeRange.Min;

                for (var k = 0; k < _isotopeList.Count; k++)
                {
                    var i = _isotopeList.SortedIndexByIntensity[k]; // internal isotope index
                    var isotopeIndex = _isotopeList[i].Index;
                    var isotopeMzLb = Ion.GetIsotopeMz(_comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex);
                    var isotopeMzUb = Ion.GetIsotopeMz(_comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex);
                    isotopeMzLb -= MzTolerance.GetToleranceAsTh(isotopeMzLb);
                    isotopeMzUb += MzTolerance.GetToleranceAsTh(isotopeMzUb);

                    if (isotopeMzLb < _run.MinMs1Mz || isotopeMzUb > _run.MaxMs1Mz) continue;

                    var st = Array.BinarySearch(_ms1PeakArr, new Peak(isotopeMzLb, 0));
                    if (st < 0) st = ~st;

                    for (var j = st; j < _ms1PeakArr.Length; j++)
                    {
                        var ms1Peak = _ms1PeakArr[j];
                        if (ms1Peak.Mz > isotopeMzUb) break;

                        var col = _ms1PeakInfoArr[j].SpectrumIndex;
                        if (k <= 3)
                        {
                            _featureMatrix[row][col][i] = Math.Max(ms1Peak.Intensity, _featureMatrix[row][col][i]);

                            if (_featureMatrix[row][col][i] > _highestPeakIntensity[row][col])
                            {
                                _highestPeakIndex[row][col] = i;
                                _highestPeakIntensity[row][col] = _featureMatrix[row][col][i];

                                if (k == 0)
                                {
                                    _mostAbuIsotopePeakIndex[row][col] = j;
                                    _accurateMass[row][col] = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                                }
                                /*else
                                {
                                    if (!(_accurateMass[row][col] > 0)) _accurateMass[row][col] = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
                                }*/
                            }
                        }
                        else
                        {
                            if (_featureMatrix[row][col][i] > 0)
                            {
                                // if existing peak matches better, then skip
                                if (_highestPeakIntensity[row][col] > 0)
                                {
                                    var expectedIntensity = (_highestPeakIntensity[row][col] * _isotopeList[i].Ratio) / _isotopeList[_highestPeakIndex[row][col]].Ratio;
                                    if (Math.Abs(ms1Peak.Intensity - expectedIntensity) < Math.Abs(_featureMatrix[row][col][i] - expectedIntensity)) // change peak
                                        _featureMatrix[row][col][i] = ms1Peak.Intensity;
                                }
                                else
                                {
                                    _featureMatrix[row][col][i] = Math.Max(ms1Peak.Intensity, _featureMatrix[row][col][i]);
                                }
                            }
                            else
                            {
                                _featureMatrix[row][col][i] = ms1Peak.Intensity;
                            }
                        }

                        if (!observedCharges[row]) observedCharges[row] = true;
                    }
                }

                if (observedCharges[row])
                {
                    for (var col = 0; col < _nScans; col++)
                    {
                        if (_highestPeakIntensity[row][col] > 0)
                        {
                            _correlationMap[row][col] = _isotopeList.GetPearsonCorrelation(_featureMatrix[row][col]);
                        }
                    }
                }
            }// end or row for-loop
            );

            var temp = new List<int>();
            for (var i = 0; i < observedCharges.Length; i++) if (observedCharges[i]) temp.Add(i);
            _chargeIndexes = temp.ToArray();
        }

        private const double CorrLowerBound = 0.5d;
        private IEnumerable<ChargeLcScanCell> GetSeedCellsByCorr()
        {
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();
            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (_correlationMap[i][j] < CorrLowerBound) continue;
                    seedCells.Add(new KeyValuePair<double, ChargeLcScanCell>(_correlationMap[i][j], new ChargeLcScanCell(i, j)));
                }
            }

            return seedCells.OrderByDescending(x => x.Key).Select(x => x.Value);
        }

        private List<ChargeLcScanCluster> FindClusters(double clusteringScoreCutoff = 0.7)
        {
            BuildFeatureMatrix(); // should be called first

            var clusters = new List<ChargeLcScanCluster>();
            var tempEnvelope = new double[_isotopeList.Count];
            int chargeNeighborGap = 2;
            const int scanNeighborGap = 1;
         
            foreach (var seedCell in GetSeedCellsByCorr())
            {
                if (_checkedOut[seedCell.Row][seedCell.Col]) continue;
                var seedScore = _correlationMap[seedCell.Row][seedCell.Col];
                var seedMass = _accurateMass[seedCell.Row][seedCell.Col];

                var newCluster = new ChargeLcScanCluster(_chargeRange.Min, _ms1ScanNums, _isotopeList);
                newCluster.AddMember(seedCell, _featureMatrix[seedCell.Row][seedCell.Col], _correlationMap[seedCell.Row][seedCell.Col]);
                newCluster.ClusteringScore2 = _isotopeList.GetBhattacharyyaDistance(_featureMatrix[seedCell.Row][seedCell.Col]);

                var neighbors = new Queue<ChargeLcScanCell>();

                neighbors.Enqueue(seedCell); // pick a seed
                _checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    
                    var charge = cell.Row + _chargeRange.Min;
                    
                    if (charge > 35) chargeNeighborGap = 8;
                    else if (charge > 25) chargeNeighborGap = 4;
                    else if (charge > 15) chargeNeighborGap = 2;

                    for (var l = Math.Max(cell.Col - scanNeighborGap, 0); l <= Math.Min(cell.Col + scanNeighborGap, _nScans - 1); l++)
                    {
                        for (var k = Math.Max(cell.Row - chargeNeighborGap, _chargeIndexes.First()); k <= Math.Min(cell.Row + chargeNeighborGap, _chargeIndexes.Last()); k++)
                        {
                            var distFromSeed = Math.Abs(seedCell.Col - l);

                            if (_checkedOut[k][l] || (_correlationMap[k][l] < CorrLowerBound && distFromSeed > 1) || (_correlationMap[k][l] < 0.2 && distFromSeed < 2)) continue;

                            if (Math.Abs(seedMass - _accurateMass[k][l]) > MzTolerance.GetToleranceAsTh(seedMass)) continue;

                            for (var t = 0; t < tempEnvelope.Length; t++) tempEnvelope[t] = newCluster.ClusteringEnvelope[t] + _featureMatrix[k][l][t];
                            var newCorr = _isotopeList.GetPearsonCorrelation(tempEnvelope);
                            var newDivergence = _isotopeList.GetBhattacharyyaDistance(tempEnvelope);

                            if ((newCluster.ClusteringScore < newCorr || newCluster.ClusteringScore2 > newDivergence) ||
                                (distFromSeed < 1 && seedScore < newCorr) || (_correlationMap[k][l] > 0.9))
                            {
                                var newMember = new ChargeLcScanCell(k, l);
                                neighbors.Enqueue(newMember);
                                newCluster.AddMember(newMember, _featureMatrix[k][l], newCorr);
                                _checkedOut[k][l] = true;
                                newCluster.ClusteringScore2 = newDivergence;
                            }
                        }
                    }
                }

                for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) _checkedOut[i][j] = true;

                if (newCluster.ClusteringScore < clusteringScoreCutoff) continue;

                var idx = _mostAbuIsotopePeakIndex[seedCell.Row][seedCell.Col];
                newCluster.SetRepresentativeMass(seedMass, _ms1PeakArr[idx].Mz, seedCell.Row + _chargeRange.Min, _ms1ScanNums[seedCell.Col]);

                clusters.Add(newCluster);
            }

            return clusters;
        }

        private List<ChargeLcScanCluster> GetFilteredClusters(List<ChargeLcScanCluster> clusters, bool filtering = true)
        {
            var options = new ParallelOptions();
            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            Parallel.ForEach(clusters, options, GetAccurateMassInfo);
            var filteredClusters = new List<ChargeLcScanCluster>();

            foreach (var cluster in clusters)
            {
                if (!cluster.Active) continue;
                if (filtering && !cluster.GoodEnough && cluster.TooBad) continue;
                filteredClusters.Add(cluster);
            }

            return filteredClusters;
        }

        private const int MinXicWindowLength = 10;
        private double CalculateXicCorrelationOverTimeBetweenIsotopes(ChargeLcScanCluster cluster)
        {
            var maxCol = cluster.MaxCol;
            var minCol = cluster.MinCol;
            var maxRow = cluster.MaxRow;
            var minRow = cluster.MinRow;
            var colLen = maxCol - minCol + 1;

            if (colLen < MinXicWindowLength)
            {
                minCol = Math.Max(minCol - (int)((MinXicWindowLength - colLen) * 0.5), 0);

                if (minCol == 0) maxCol = minCol + MinXicWindowLength - 1;
                else
                {
                    maxCol = Math.Min(maxCol + (int)((MinXicWindowLength - colLen) * 0.5), _nScans - 1);
                    if (maxCol == _nScans - 1) minCol = maxCol - MinXicWindowLength + 1;
                }
                colLen = maxCol - minCol + 1;
            }

            var n = Math.Min(_isotopeList.Count, 3);
            var xicProfile = new double[n][];
            for (var i = 0; i < xicProfile.Length; i++) xicProfile[i] = new double[colLen];

            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    for (var k = 0; k < n; k++)
                    {
                        var isoIndex = _isotopeList.SortedIndexByIntensity[k];
                        xicProfile[k][col - minCol] += _featureMatrix[row][col][isoIndex];
                    }
                }
            }

            // smoothing
            //for (var k = 0; k < _topEnvelopes.Length; k++) xicProfile[k] = _smoother.Smooth(xicProfile[k]);
            if (n < 3) return FitScoreCalculator.GetPearsonCorrelation(xicProfile[0], xicProfile[1]);

            double ret = Math.Max(FitScoreCalculator.GetPearsonCorrelation(xicProfile[0], xicProfile[1]),
                                FitScoreCalculator.GetPearsonCorrelation(xicProfile[0], xicProfile[2]));
            ret = Math.Max(ret, FitScoreCalculator.GetPearsonCorrelation(xicProfile[2], xicProfile[1]));

            return ret;
        }

        private void GetAccurateMassInfo(ChargeLcScanCluster cluster)
        {
            var bestEnvelopeCorrelation = -1.0d;
            var bestRankSumScore = 1.0d;
            var bestPoissonScore = 1.0d;
            var bestBhattacharyyaDistance = 999999d;
            var bestKullbackLeiblerDivergence = 999999d;

            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;

            var memberCorr = new List<double>();
            var memberRankSum = new List<double>();
            var memberPoisson = new List<double>();
            var memberBc = new List<double>();
            //var memberkl = new List<double>();

            var summedN = 0;
            var summedK = 0;
            var summedN1 = 0;
            var summedK1 = 0;

            var memberN = new List<int>();
            var memberK = new List<int>();
            var memberN1 = new List<int>();
            var memberK1 = new List<int>();
            var memberObsPeak = new List<Peak[]>();
            
            var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);

            cluster.ClearMember();

            for (var col = minCol; col <= maxCol; col++)
            {
                for (var row = minRow; row <= maxRow; row++)
                {
                    var mass = _accurateMass[row][col];

                    if (!(mass > 0)) continue;
                    if (Math.Abs(cluster.RepresentativeMass - mass) > massTol) continue;

                    var charge = row + _chargeRange.Min;
                    var envelope = new double[_isotopeList.Count];
                    var spectrum = _cachedSpectrum[col];
                    var mostAbuPeakIndex = _mostAbuIsotopePeakIndex[row][col];
                    int nPeaks;
                    Peak[] observedPeaks;
                    int[] observedPeakRanks;
                    int nCheckedIsotopes;
                    int nObservedIsotopes;
                    double binMinMz;
                    double binMaxMz;

                    var mostAbuPeakIndexInSpectrum = _ms1PeakInfoArr[mostAbuPeakIndex].PeakIndexInSpectrum;
                    
                    spectrum.GetRankSum(_isotopeList, charge, mass, mostAbuPeakIndexInSpectrum, out nPeaks, out observedPeaks, out observedPeakRanks, out nCheckedIsotopes, out nObservedIsotopes, out binMinMz, out binMaxMz);

                    var n = (int)Math.Round((binMaxMz - binMinMz) / (0.5 * (MzTolerance.GetToleranceAsTh(binMaxMz) + MzTolerance.GetToleranceAsTh(binMinMz))));
                    var k = nCheckedIsotopes; // # of theretical isotope ions of the mass within the local window
                    var n1 = nPeaks; // # of detected ions within the local window
                    var k1 = nObservedIsotopes; // # of matched ions for generating isotope envelope profile

                    if (nObservedIsotopes < 3) continue;

                    for (var i = 0; i < _isotopeList.Count; i++)
                    {
                        if (observedPeaks[i] == null) continue;
                        envelope[i] = observedPeaks[i].Intensity;
                    }

                    var envCorr = _isotopeList.GetPearsonCorrelation(envelope);

                    memberN.Add(n);
                    memberK.Add(k);
                    memberN1.Add(n1);
                    memberK1.Add(k1);
                    memberObsPeak.Add(observedPeaks);

                    var ranksumScore = CalculateRankSumScore(nPeaks, observedPeakRanks);
                    var poissonScore = CalculatePoissonScore(n, k, n1, k1);

                    var bcDistance = _isotopeList.GetBhattacharyyaDistance(envelope);
                    var klDivergence = _isotopeList.GetKullbackLeiblerDivergence(envelope);

                    if (bestPoissonScore < poissonScore) bestPoissonScore = poissonScore;
                    if (bestRankSumScore < ranksumScore) bestRankSumScore = ranksumScore;
                    if (bestEnvelopeCorrelation < envCorr) bestEnvelopeCorrelation = envCorr;
                    if (bcDistance < bestBhattacharyyaDistance) bestBhattacharyyaDistance = bcDistance;
                    if (klDivergence < bestKullbackLeiblerDivergence) bestKullbackLeiblerDivergence = klDivergence;

                    cluster.AddMember(new ChargeLcScanCell(row, col), envelope);

                    memberRankSum.Add(ranksumScore);
                    memberCorr.Add(envCorr);
                    memberPoisson.Add(poissonScore);
                    memberBc.Add(bcDistance);
                }
            }

            // pre-filtering condition
            if (cluster.Members.Count < 1)
            {
                cluster.Active = false;
                return;
            }

            cluster.SetScore(ChargeLcScanScore.EnvelopeCorrelation, bestEnvelopeCorrelation);
            cluster.SetScore(ChargeLcScanScore.RankSum, bestRankSumScore);
            cluster.SetScore(ChargeLcScanScore.Poisson, bestPoissonScore);
            cluster.SetScore(ChargeLcScanScore.BhattacharyyaDistance, bestBhattacharyyaDistance);
            cluster.SetScore(ChargeLcScanScore.KullbackLeiblerDivergence, bestKullbackLeiblerDivergence);

            double summedBc;
            double summedCorr;
            double summedKl;
            List<int> selectedMembers;
            var summedEnvelope = GetSummedEnvelope(cluster, memberCorr, memberBc, out summedCorr, out summedBc, out summedKl, out selectedMembers);

            cluster.SummedEnvelope = summedEnvelope;
            cluster.SetScore(ChargeLcScanScore.EnvelopeCorrelationSummed, summedCorr);
            cluster.SetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed, summedBc);
            cluster.SetScore(ChargeLcScanScore.KullbackLeiblerDivergenceSummed, summedKl);

            var avgRankSum = 0d;
            var avgPoission = 0d;
            var totalMzErrorPpm = 0d;
            var nPair = 0;
            var bestErrorPpm = 999999d;
            foreach (var i in selectedMembers)
            {
                avgRankSum += memberRankSum[i];
                avgPoission += memberPoisson[i];
                summedN += memberN[i];
                summedK += memberK[i];
                summedN1 += memberN1[i];
                summedK1 += memberK1[i];

                var charge = cluster.Members[i].Row + _chargeRange.Min;
                var peaks = memberObsPeak[i];
                var n = 0;
                var mzErrorPpm = 0d;
                for (var j = 0; j < _isotopeList.Count; j++)
                {
                    if (peaks[j] == null) continue;
                    var theoreticalMz = Ion.GetIsotopeMz(cluster.RepresentativeMass, charge, _isotopeList[j].Index);
                    mzErrorPpm += (Math.Abs(peaks[j].Mz - theoreticalMz) * 1e6) / theoreticalMz;
                    n++;
                }
                totalMzErrorPpm += mzErrorPpm;
                nPair += n;
                mzErrorPpm /= n;

                if (bestErrorPpm > mzErrorPpm) bestErrorPpm = mzErrorPpm;
            }

            totalMzErrorPpm /= nPair;
            avgRankSum /= selectedMembers.Count;
            avgPoission /= selectedMembers.Count;
            cluster.SetScore(ChargeLcScanScore.MzError, bestErrorPpm);
            cluster.SetScore(ChargeLcScanScore.MzErrorSummed, totalMzErrorPpm);
            cluster.SetScore(ChargeLcScanScore.RankSumMedian, avgRankSum);
            cluster.SetScore(ChargeLcScanScore.PoissonMedian, avgPoission);
            cluster.SetScore(ChargeLcScanScore.PoissonSummed, CalculatePoissonScore(summedN, summedK, summedN1, summedK1));

            var xicCorr = CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
            cluster.SetScore(ChargeLcScanScore.XicCorrelation, xicCorr);

            cluster.SetBoundary(selectedMembers);
            cluster.Probability = cluster.GetProbabilityByLogisticRegression();
        }

        private double[] GetSummedEnvelope(ChargeLcScanCluster cluster, List<double> cellCorr, List<double> cellBc, out double summedCorr, out double summedBc, out double summedKl, out List<int> selectedMembers)
        {
            var summedEnvelop = new double[_isotopeList.Count];
            var tempEnvelop = new double[_isotopeList.Count];
            var seleted = new bool[cellCorr.Count];

            summedBc = 99999d;
            var index2 = Enumerable.Range(0, cellBc.Count).ToArray();
            var arr2 = cellBc.ToArray();
            Array.Sort(arr2, index2);
            foreach (var j in index2)
            {
                Array.Clear(tempEnvelop, 0, tempEnvelop.Length);
                for (var k = 0; k < _isotopeList.Count; k++) tempEnvelop[k] = summedEnvelop[k] + cluster.MemberEnvelope[j][k];

                var tempBc = _isotopeList.GetBhattacharyyaDistance(tempEnvelop);

                if (tempBc > summedBc) continue;
                summedBc = tempBc;

                seleted[j] = true;
                Array.Copy(tempEnvelop, summedEnvelop, tempEnvelop.Length);
            }

            summedKl = _isotopeList.GetKullbackLeiblerDivergence(summedEnvelop);

            var index = Enumerable.Range(0, cellCorr.Count).ToArray();
            var arr = cellCorr.ToArray();
            Array.Sort(arr, index);
            Array.Reverse(index);
            summedCorr = 0d;
            Array.Clear(summedEnvelop, 0, summedEnvelop.Length);
            foreach (var j in index)
            {
                Array.Clear(tempEnvelop, 0, tempEnvelop.Length);
                
                for (var k = 0; k < _isotopeList.Count; k++) tempEnvelop[k] = summedEnvelop[k] + cluster.MemberEnvelope[j][k];
                var tempCorr = _isotopeList.GetPearsonCorrelation(tempEnvelop);

                if (tempCorr < summedCorr)
                {
                     if (_isotopeList.GetPearsonCorrelation(cluster.MemberEnvelope[j]) < 0.9) continue;
                }
                
                summedCorr = tempCorr;
                
                seleted[j] = true;
                Array.Copy(tempEnvelop, summedEnvelop, tempEnvelop.Length);
            }

            selectedMembers = new List<int>();
            for (var i = 0; i < seleted.Length; i++)
                if(seleted[i]) selectedMembers.Add(i);
            
            cluster.EnvelopeCount = selectedMembers.Count;

            return summedEnvelop;
        }

        private static double CalculatePoissonScore(int n, int k, int n1, int k1)
        {
            var lambda = ((double)n1 / (double)n) * k;
            var pvalue = 1 - Poisson.CDF(lambda, k1);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        private double CalculateRankSumScore(int nPeaks, int[] isotopePeakRanks)
        {
            var n = 0;
            var r = 0d;
            for (var k = 0; k < _isotopeList.Count; k++)
            {
                var internalIndex = _isotopeList.SortedIndexByIntensity[k];
                if (_isotopeList[internalIndex].Ratio < 0.7) break;

                if (isotopePeakRanks[internalIndex] > 0)
                {
                    n++;
                    r += isotopePeakRanks[internalIndex];
                }
            }

            if (nPeaks == n) return 0d;
            var pvalue = FitScoreCalculator.GetRankSumPvalue(nPeaks, n, r);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        public double[][] GetCorrelationMap()
        {
            return _correlationMap;
        }

        private int _minSearchMassBin;
        private int _maxSearchMassBin;
        List<ChargeLcScanCluster>[] _massBinToClusterMap;
        private void ExtractFeatures(int binNum, bool massCollapse, double probabilityThreshold)
        {
            var monoIsotopicMass = _comparer.GetMzAverage(binNum);
            var neighborMassBins = new List<int>();
            if (binNum > _minSearchMassBin) neighborMassBins.Add(binNum - 1);
            if (binNum < _maxSearchMassBin) neighborMassBins.Add(binNum + 1);

            if (massCollapse)
            {
                for (var i = -2; i <= 2; i++)
                {
                    if (i == 0) continue;
                    var neighborIsotopebinNum = _comparer.GetBinNumber(monoIsotopicMass + i);
                    if (neighborIsotopebinNum >= _minSearchMassBin && neighborIsotopebinNum <= _maxSearchMassBin) neighborMassBins.Add(neighborIsotopebinNum);
                }
            }

            foreach (var cluster in GetProbableChargeScanClusters(monoIsotopicMass))
            {
                if (cluster.Probability < probabilityThreshold && !cluster.GoodEnough) continue;

                // Before adding it, Check if there are existing neighbor that can be merged
                // search mass range = +-2 [Da]
                
                var foundNeighbor = false;
                //var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);
                foreach (var neighborBin in neighborMassBins)
                {
                    foreach (var neighborCluster in _massBinToClusterMap[neighborBin - _minSearchMassBin])
                    {
                        foundNeighbor = ChargeLcScanCluster.Merge(cluster, neighborCluster);
                    }

                    if (foundNeighbor) break;
                }
                _massBinToClusterMap[binNum - _minSearchMassBin].Add(cluster);
            }
        }

        private int _nOutFatures;
        private void FlushOutput(int startBin, int endBin, StreamWriter tsvWriter, double probabilityThreshold, bool scoreReport, StreamWriter csvWriter)
        {
            for (var binNum = startBin; binNum <= endBin; binNum++)
            {
                foreach (var cluster in _massBinToClusterMap[binNum - _minSearchMassBin])
                {
                    if (!cluster.Active) continue;
                    if (cluster.GoodEnough || cluster.Probability > probabilityThreshold)
                    {
                        tsvWriter.WriteLine(cluster.GetString(scoreReport));
                        _nOutFatures++;

                        if (csvWriter != null)
                        {
                            //csvWriter.WriteLine("scan_num,charge,monoisotopic_mw,mz,abundance");
                            for (var j = 0; j < cluster.Members.Count; j++)
                            {
                                var m = cluster.Members[j];
                                var idx = _mostAbuIsotopePeakIndex[m.Row][m.Col];
                                var mz = _ms1PeakArr[idx].Mz;
                                var abundance = cluster.MemberEnvelope[j].Sum();
                                csvWriter.WriteLine(string.Format("{0},{1},{2},{3},{4},{5}", _ms1ScanNums[m.Col], m.Row+_chargeRange.Min, _accurateMass[m.Row][m.Col], mz, abundance, cluster.Probability));
                            }
                        }
                    }
                }
                _massBinToClusterMap[binNum - _minSearchMassBin].Clear();
            }
        }

    }

    public class IsotopeList : List<Isotope>
    {
        public double MonoIsotopeMass { get; private set; }
        private readonly double[] _envelope;
        private readonly double[] _envelopePdf;
        public int[] SortedIndexByIntensity { get; private set; }

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

            _envelope = this.Select(iso => iso.Ratio).ToArray();
            SortedIndexByIntensity = new int[Count];
            _envelopePdf = new double[Count];
            var s = _envelope.Sum();

            for (var i = 0; i < Count; i++)
            {
                var rankingIndex = isotopeRankings[this[i].Index] - 1;
                SortedIndexByIntensity[rankingIndex] = i;

                _envelopePdf[i] = _envelope[i] / s;
            }
        }

        public Isotope GetIsotopeRankedAt(int ranking)
        {
            return this[SortedIndexByIntensity[ranking - 1]];
        }

        public double GetPearsonCorrelation(double[] observedIsotopeEnvelop)
        {
            var m1 = 0.0;
            var m2 = 0.0;
            var nMissedPeaks = 0;

            for (var i = 0; i < Count; i++)
            {
                m1 += _envelope[i];
                if (observedIsotopeEnvelop[i] > 0) m2 += observedIsotopeEnvelop[i];
                else nMissedPeaks++;
            }
            
            /*
            if (MonoIsotopeMass < 13000)
            {
                if (nMissedPeaks > Count * 0.7) return 0.0d;
            }
            else
            {
                if (nMissedPeaks > Count * 0.7) return 0.0d;
            }
            */

            m1 /= Count;
            m2 /= Count;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < Count; i++)
            {
                var d1 = _envelope[i] - m1;
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
            return FitScoreCalculator.GetBhattacharyyaDistance(_envelope, observedIsotopeEnvelop, _envelope.Length);
        }

        public double GetJensenShannonDivergence(double[] observedIsotopeEnvelop)
        {
            var n = 0;
            var s = 0d;
            var ret = 0d;
            for (var i = 0; i < Count; i++)
            {
                s += observedIsotopeEnvelop[i];
                if (observedIsotopeEnvelop[i] > 0) n++;
            }

            if (n != Count)
            {
                const double eps = 1e-6;
                var qc = eps / n;

                for (var i = 0; i < Count; i++)
                {
                    var observedProb = observedIsotopeEnvelop[i] > 0 ? observedIsotopeEnvelop[i] / s - qc : eps;
                    var prob = 0.5 * (_envelopePdf[i] + observedProb);
                    ret += (observedProb * (Math.Log(observedProb, 2) - Math.Log(prob, 2))) * 0.5;
                    ret += (_envelopePdf[i] * (Math.Log(_envelopePdf[i], 2) - Math.Log(prob, 2))) * 0.5;
                }
            }
            else
            {
                for (var i = 0; i < Count; i++)
                {
                    var observedProb = observedIsotopeEnvelop[i] / s;
                    var prob = 0.5 * (_envelopePdf[i] + observedProb);
                    ret += (observedProb * (Math.Log(observedProb, 2) - Math.Log(prob, 2))) * 0.5;
                    ret += (_envelopePdf[i] * (Math.Log(_envelopePdf[i], 2) - Math.Log(prob, 2))) * 0.5;
                }
            }

            return ret;
        }

        public double GetKullbackLeiblerDivergence(double[] observedIsotopeEnvelop, bool smoothing = true)
        {
            const double eps = 1e-6;
            var n = 0;
            var s = 0d;
            var ret = 0d;

            for (var i = 0; i < Count; i++)
            {
                s += observedIsotopeEnvelop[i];
                if (observedIsotopeEnvelop[i] > 0) n++;
            }

            if (n != Count)
            {
                var qc = eps / n;
                for (var i = 0; i < Count; i++)
                {
                    var observedProb = observedIsotopeEnvelop[i] > 0 ? observedIsotopeEnvelop[i] / s - qc : eps;
                    ret += _envelopePdf[i] * (Math.Log(_envelopePdf[i], 2) - Math.Log(observedProb, 2));
                }
            }
            else
            {
                for (var i = 0; i < Count; i++)
                {
                    var observedProb = observedIsotopeEnvelop[i] / s;
                    ret += _envelopePdf[i] * (Math.Log(_envelopePdf[i], 2) - Math.Log(observedProb, 2));
                }

            }
            return ret;
        }

        public Isotope GetMostAbundantIsotope()
        {
            return GetIsotopeRankedAt(1);
        }
    }

    class Ms1PeakInfo
    {
        internal int SpectrumIndex { get; private set; }
        internal int PeakIndexInSpectrum { get; private set; }

        internal Ms1PeakInfo(int specIndex, int peakIndex)
        {
            SpectrumIndex = specIndex;
            PeakIndexInSpectrum = peakIndex;
        }
    }

    class ChargeLcScanSpectrum : Spectrum
    {
        internal ChargeLcScanSpectrum(Spectrum spec, int numBits4WinSize = 18)
            : base(spec.Peaks, spec.ScanNum)
        {
            WindowComparer = new MzComparerWithBinning(numBits4WinSize);
            _minBinNum = WindowComparer.GetBinNumber(MinMz);
            _maxBinNum = WindowComparer.GetBinNumber(MaxMz);
            _peakIndex = new IntRange[3][];
            _peakRanking = new int[3][][];
            for (var i = 0; i < 3; i++)
            {
                _peakIndex[i] = new IntRange[_maxBinNum - _minBinNum + 1];
                _peakRanking[i] = new int[_maxBinNum - _minBinNum + 1][];
            }

            for (var binNum = _minBinNum; binNum <= _maxBinNum; binNum++)
            {
                var minMz = WindowComparer.GetMzAverage(binNum - 1);
                var maxMz = WindowComparer.GetMzAverage(binNum + 1);

                var startIndex = Array.BinarySearch(Peaks, new Peak(minMz - float.Epsilon, 0));
                if (startIndex < 0) startIndex = ~startIndex;

                // given a window, scan spectrum from three different start position s.t. accurate coverage
                var intensities0 = new List<double>();
                var minMz0 = WindowComparer.GetMzStart(binNum);
                var maxMz0 = WindowComparer.GetMzEnd(binNum);

                var intensities1 = new List<double>();
                var minMz1 = WindowComparer.GetMzAverage(binNum - 1);
                var maxMz1 = WindowComparer.GetMzAverage(binNum);

                var intensities2 = new List<double>();
                var minMz2 = maxMz1;
                var maxMz2 = WindowComparer.GetMzAverage(binNum + 1);

                var start0 = Peaks.Length;
                var start1 = Peaks.Length;
                var start2 = Peaks.Length;
                var end0 = -1;
                var end1 = -1;
                var end2 = -1;

                var i = startIndex;
                while (i < Peaks.Length)
                {
                    if (Peaks[i].Mz > maxMz) break;

                    if (Peaks[i].Mz >= minMz0 && Peaks[i].Mz < maxMz0)
                    {
                        intensities0.Add(Peaks[i].Intensity);
                        if (i < start0) start0 = i;
                        if (i > end0) end0 = i;
                    }
                    if (Peaks[i].Mz >= minMz1 && Peaks[i].Mz < maxMz1)
                    {
                        intensities1.Add(Peaks[i].Intensity);
                        if (i < start1) start1 = i;
                        if (i > end1) end1 = i;
                    }
                    if (Peaks[i].Mz >= minMz2 && Peaks[i].Mz < maxMz2)
                    {
                        intensities2.Add(Peaks[i].Intensity);
                        if (i < start2) start2 = i;
                        if (i > end2) end2 = i;
                    }
                    ++i;
                }

                if (end0 > -1)
                {
                    _peakIndex[0][binNum - _minBinNum] = new IntRange(start0, end0);
                    _peakRanking[0][binNum - _minBinNum] = ArrayUtil.GetRankings(intensities0);
                }
                if (end1 > -1)
                {
                    _peakIndex[1][binNum - _minBinNum] = new IntRange(start1, end1);
                    _peakRanking[1][binNum - _minBinNum] = ArrayUtil.GetRankings(intensities1);
                }

                if (end2 > -1)
                {
                    _peakIndex[2][binNum - _minBinNum] = new IntRange(start2, end2);
                    _peakRanking[2][binNum - _minBinNum] = ArrayUtil.GetRankings(intensities2);
                }
            }
        }

        internal MzComparerWithBinning WindowComparer { private set; get; }

        internal int GetRanks(int binShift, int binIndex, int globalIndex)
        {
            var ranking = 0;
            if (globalIndex >= _peakIndex[binShift][binIndex].Min && globalIndex <= _peakIndex[binShift][binIndex].Max)
                ranking = _peakRanking[binShift][binIndex][globalIndex - _peakIndex[binShift][binIndex].Min];

            return ranking;
        }

        private void GetWindow(double mz, out int binShift, out int binIndex, out double binMzStart, out double binMzEnd)
        {
            var binNum = WindowComparer.GetBinNumber(mz);
            binIndex = binNum - _minBinNum;
            binShift = 0;

            var d0 = Math.Abs(WindowComparer.GetMzAverage(binNum) - mz);
            var d1 = Math.Abs(WindowComparer.GetMzStart(binNum) - mz);
            var d2 = Math.Abs(WindowComparer.GetMzEnd(binNum) - mz);

            binMzStart = WindowComparer.GetMzStart(binNum);
            binMzEnd = WindowComparer.GetMzEnd(binNum);

            if (d1 < d2 && d1 < d0)
            {
                binShift = 1;
                binMzStart = WindowComparer.GetMzAverage(binNum - 1);
                binMzEnd = WindowComparer.GetMzAverage(binNum);
            }
            else if (d2 < d1 && d2 < d0)
            {
                binShift = 2;
                binMzStart = WindowComparer.GetMzAverage(binNum);
                binMzEnd = WindowComparer.GetMzAverage(binNum + 1);
            }

        }

        internal void GetRankSum(IsotopeList isotopeList, int charge, double monoIsotopeMass, int mostAbundantIsotopePeakIndex,
                                out int nPeaks, out Peak[] observedPeaks, out int[] observedPeakRanks,
                                out int nCheckedIsotopes, out int nObservedIsotopes, out double binMzStart, out double binMzEnd)
        {
            var tolerance = ChargeLcScanMatrix.MzTolerance;
            var mostAbundantIsotopeMz = Peaks[mostAbundantIsotopePeakIndex].Mz;

            int binShift;
            int binIndex;
            GetWindow(mostAbundantIsotopeMz, out binShift, out binIndex, out binMzStart, out binMzEnd);

            if (_peakIndex[binShift][binIndex] == null) throw new ArgumentOutOfRangeException();

            nPeaks = _peakIndex[binShift][binIndex].Max - _peakIndex[binShift][binIndex].Min + 1;

            var mostAbuInternalIndex = isotopeList.SortedIndexByIntensity[0];
            observedPeaks = new Peak[isotopeList.Count];
            observedPeakRanks = new int[isotopeList.Count];
            observedPeaks[mostAbuInternalIndex] = Peaks[mostAbundantIsotopePeakIndex];
            observedPeakRanks[mostAbuInternalIndex] = GetRanks(binShift, binIndex, mostAbundantIsotopePeakIndex);

            // # of checked isotopes within local window
            nCheckedIsotopes = 1;
            nObservedIsotopes = 1;

            // go down
            var peakIndex = mostAbundantIsotopePeakIndex - 1;
            for (var isotopeInternalIndex = isotopeList.SortedIndexByIntensity[0] - 1; isotopeInternalIndex >= 0; isotopeInternalIndex--)
            {
                var isotopeIndex = isotopeList[isotopeInternalIndex].Index; //selectedIsotopeIndex[isotopeInternalIndex];
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;

                if (binMzStart <= isotopeMz && isotopeMz <= binMzEnd) nCheckedIsotopes++;
                else break;

                for (var i = peakIndex; i >= 0; i--)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        var peak = Peaks[i];
                        if (observedPeaks[isotopeInternalIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeInternalIndex].Intensity)
                        {
                            observedPeaks[isotopeInternalIndex] = peak;
                            observedPeakRanks[isotopeInternalIndex] = GetRanks(binShift, binIndex, i);
                            if (observedPeakRanks[isotopeInternalIndex] > 0) nObservedIsotopes++;
                        }
                    }
                }
            }

            // go up
            peakIndex = mostAbundantIsotopePeakIndex + 1;
            for (var isotopeInternalIndex = mostAbuInternalIndex + 1; isotopeInternalIndex < isotopeList.Count; isotopeInternalIndex++)
            {
                var isotopeIndex = isotopeList[isotopeInternalIndex].Index;
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;

                if (binMzStart <= isotopeMz && isotopeMz <= binMzEnd) nCheckedIsotopes++;
                else break;

                for (var i = peakIndex; i < Peaks.Length; i++)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz >= minMz)    // find match, move to next isotope
                    {
                        var peak = Peaks[i];
                        if (observedPeaks[isotopeInternalIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeInternalIndex].Intensity)
                        {
                            observedPeaks[isotopeInternalIndex] = peak;
                            observedPeakRanks[isotopeInternalIndex] = GetRanks(binShift, binIndex, i);
                            if (observedPeakRanks[isotopeInternalIndex] > 0) nObservedIsotopes++;
                        }
                    }
                }
            }
        }

        internal double MinMz { get { return Peaks[0].Mz; } }
        internal double MaxMz { get { return Peaks[Peaks.Length - 1].Mz; } }

        private readonly IntRange[][] _peakIndex;
        private readonly int[][][] _peakRanking;

        private readonly int _minBinNum;
        private readonly int _maxBinNum;
    }
}
