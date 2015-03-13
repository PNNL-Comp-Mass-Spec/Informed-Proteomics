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
    public class Ms1FeatureMatrix : IMs1FeatureExtract
    {
        public Ms1FeatureMatrix(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0, int numBits = 27)
        {
            _run = run;
            _minScanCharge  = minScanCharge;
            _maxScanCharge  = maxScanCharge;
            _maxThreadCount = maxThreadCount;
            Comparer       = new MzComparerWithBinning(numBits);
            _ms1ScanNums    = run.GetMs1ScanVector();
            _nScans         = _ms1ScanNums.Length;
            _ms1PeakList    = new List<Ms1Peak>();
            _spectrums      = new List<Ms1Spectrum>();

            Ms1FeatureCluster.HasMs2Spectra = (run.MaxMsLevel > 1);
            //Console.WriteLine("Number of ms1 spectra = {0}", _ms1ScanNums.Length);
            for (var i = 0; i < _ms1ScanNums.Length; i++)
            {
                //var spec = _run.GetSpectrum(_ms1ScanNums[i]);
                //var ms1Spec = new Ms1Spectrum(spec.ScanNum, spec.Peaks, i);
                var ms1Spec = run.GetMs1Spectrum(_ms1ScanNums[i]);
                _spectrums.Add(ms1Spec);
                _ms1PeakList.AddRange((Ms1Peak[])ms1Spec.Peaks);
                //Console.WriteLine("{0}th ms1 spectrum is loaded; # of peaks = {1}", i, ms1Spec.Peaks.Length);
            }
            _ms1PeakList.Sort();

            _correlationMap     = new double[MaxChargeLength][];
            _distanceMap        = new double[MaxChargeLength][];
            _accurateMass       = new double[MaxChargeLength][];
            _featureMatrix      = new Ms1Peak[MaxChargeLength][][];
            _checkedOut         = new bool[MaxChargeLength][];

            for (var i = 0; i < MaxChargeLength; i++)
            {
                _checkedOut[i]          = new bool[_nScans];
                _correlationMap[i]      = new double[_nScans];
                _featureMatrix[i]       = new Ms1Peak[_nScans][];
                _distanceMap[i]         = new double[_nScans];
                _accurateMass[i]        = new double[_nScans];

                for (var j = 0; j < _nScans; j++)
                {
                    _featureMatrix[i][j] = new Ms1Peak[MaxEnvelopeLength];
                }
            }
        }

        public string GetFeatureFile(string rawFilePath, double minSearchMass = 3000, double maxSearchMass = 50000)
        {
            var outTsvFilePath = GetFeatureFilePath(rawFilePath);
            if (File.Exists(outTsvFilePath)) return outTsvFilePath;
            return GenerateFeatureFile(rawFilePath, minSearchMass, maxSearchMass, false, true);
        }

        private int _minSearchMassBin;
        private int _maxSearchMassBin;
        
        public string GenerateFeatureFile(string rawFilePath, double minMass = 3000, double maxMass = 50000, bool scoreReport = false, bool csvOutput = false)
        {
            var j = rawFilePath.LastIndexOf('.');
            var outPath = rawFilePath.Substring(0, j);

            var outTsvFilePath = GetFeatureFilePath(rawFilePath);
            var outCsvFilePath = string.Format("{0}_ms1ft.csv", outPath);
            var container = new Ms1FeatureContainer(_spectrums);

            if (File.Exists(outTsvFilePath)) return outTsvFilePath;

            _minSearchMassBin = Comparer.GetBinNumber(minMass);
            _maxSearchMassBin = Comparer.GetBinNumber(maxMass);

            double totalMassBin = _maxSearchMassBin - _minSearchMassBin + 1;

            Console.WriteLine("Start MS1 feature extracting...");
            //Console.WriteLine("Mass Range: {0} - {1}", minMass, maxMass);
            //Console.WriteLine("Charge Range: {0} - {1}", _minScanCharge, _maxScanCharge);
            Console.WriteLine("Output File: {0}", outTsvFilePath);
            if (csvOutput) Console.WriteLine("Csv Output File\t{0}", outCsvFilePath);

            var stopwatch = Stopwatch.StartNew();
            for (var binNum = _minSearchMassBin; binNum <= _maxSearchMassBin; binNum++)
            {
                var clusters = GetProbableChargeScanClusters(binNum);
                container.Add(clusters);

                if (binNum > _minSearchMassBin && (binNum - _minSearchMassBin) % 3000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - _minSearchMassBin;
                    var processedPercentage = ((double)processedBins / totalMassBin) * 100;
                    Console.WriteLine("Processing {0:0.00}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}", processedPercentage, Comparer.GetMzEnd(binNum), elapsed, container.NumberOfFeatures);
                }
            }
            Console.WriteLine("Completed MS1 feature finding; elapsed time = {0:0.000} sec; # of features = {1}", (stopwatch.ElapsedMilliseconds) / 1000.0d, container.NumberOfFeatures);

            StreamWriter csvWriter = null;
            if (csvOutput)
            {
                csvWriter = new StreamWriter(outCsvFilePath);
                csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
            }

            Console.WriteLine("Selecting mutually independent features from feature network graph");
            //Console.WriteLine("Generating feature network graphs from {0} features", container.NumberOfFeatures);
            var connectedFeatures = container.GetAllConnectedFeatures();
            var nCc = connectedFeatures.Count;

            //Console.WriteLine("Number of connected components = {0};  elapsed time = {1:0.000} sec", connectedFeatures.Count, (stopwatch.ElapsedMilliseconds) / 1000.0d);
            /*
            var tmpTsvFilePath = string.Format("{0}_ms1ft.tmp", outPath);
            var tmpTsvWriter = new StreamWriter(tmpTsvFilePath);
            tmpTsvWriter.WriteLine(Ms1FeatureCluster.GetHeaderString(scoreReport));
            var clusterId = 0;
            foreach (var featureSet in connectedFeatures)
            {
                clusterId++;
                foreach (var feature in featureSet)
                {
                    tmpTsvWriter.WriteLine("{0}\t{1}", clusterId, feature.GetString(scoreReport));                    
                }
            }
            tmpTsvWriter.Close();
            */
            var filteredFeatures = container.GetFilteredFeatures(connectedFeatures);
            Console.WriteLine("# of selected features = {0}; # of connected components = {1}; elapsed time = {2:0.000} sec", filteredFeatures.Count, nCc, (stopwatch.ElapsedMilliseconds) / 1000.0d);
            stopwatch.Stop();

            var featureId = 0;
            var tsvWriter = new StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(Ms1FeatureCluster.GetHeaderString(scoreReport));
            foreach (var cluster in filteredFeatures)
            {
                featureId++;
                tsvWriter.WriteLine("{0}\t{1}", featureId, cluster.GetString(scoreReport));
                
                if (csvWriter != null)
                {
                    foreach(var envelope in cluster.Envelopes)
                    {
                        var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                        var mostAbuPeak = envelope.Peaks[mostAbuIsotopeInternalIndex];
                        var charge = envelope.Row + _chargeRange.Min;
                        var mass = Ion.GetMonoIsotopicMass(mostAbuPeak.Mz, charge, cluster.IsotopeList[mostAbuIsotopeInternalIndex].Index);
                        var abundance = envelope.Abundance;
                        csvWriter.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6}", _ms1ScanNums[envelope.Col], charge, abundance, mostAbuPeak.Mz, 1.0 - cluster.Probability, mass, featureId));
                    }
                }
            }
            tsvWriter.Close();

            Console.WriteLine("Completed feature extraction");
            Console.WriteLine("ProMex output: {0}", outTsvFilePath);
            if (csvOutput)
            {
                csvWriter.Close();
                Console.WriteLine("ProMex output in ICR2LS format: {0}", outCsvFilePath);
            }

            //File.Move(tmpTsvFilePath, outTsvFilePath);
            return outTsvFilePath;
        }

        public IList<Ms1FeatureCluster> GetProbableChargeScanClusters(int queryMassBinNum)
        {
            return GetProbableChargeScanClusters(Comparer.GetMzAverage(queryMassBinNum));
        }
        
        public IList<Ms1FeatureCluster> GetProbableChargeScanClusters(double queryMass)
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
            for (var charge = _minScanCharge; charge <= _maxScanCharge; charge++)
            {
                var mostAbuMz = Ion.GetIsotopeMz(cluster.RepresentativeMass, charge, _isotopeList.GetMostAbundantIsotope().Index);
                var ms2ScanNums = _run.GetFragmentationSpectraScanNums(mostAbuMz);
                ms2ScanNumSet.AddRange(ms2ScanNums.Where(sn => sn > cluster.MinScanNum && sn < cluster.MaxScanNum));
            }

            return ms2ScanNumSet.Distinct();
        }

        private static LcMsRun _run;
        private static int _minScanCharge;
        private static int _maxScanCharge;
        private static int _nScans;
        private static int[] _ms1ScanNums;
        private static int _maxThreadCount;
        public readonly MzComparerWithBinning Comparer;

        //internal readonly static Tolerance MzTolerance = new Tolerance(10);
        internal readonly static Tolerance MzTolerance = new Tolerance(5);
        
        private static List<Ms1Spectrum> _spectrums;
        private static List<Ms1Peak> _ms1PeakList;
        private const int MaxEnvelopeLength = 30;
        private const int MaxChargeLength = 40;

        private readonly double[][] _correlationMap;
        private readonly double[][] _distanceMap;
        private readonly double[][] _accurateMass;
        private readonly bool[][] _checkedOut;
        private readonly Ms1Peak[][][] _featureMatrix;

        private IsotopeList _isotopeList;
        private int[] _chargeIndexes;
        private IntRange _chargeRange;
        private double _queryMass;

        public static string GetFeatureFilePath(string rawFilePath)
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
            var queryMassBinNum     = Comparer.GetBinNumber(_queryMass);
            var options             = new ParallelOptions();
            var observedCharges     = new bool[_chargeRange.Length];
            var rows                = Enumerable.Range(0, _chargeRange.Length);

            _isotopeList = new IsotopeList(Comparer.GetMzAverage(queryMassBinNum), MaxEnvelopeLength);

            if (_maxThreadCount > 0) options.MaxDegreeOfParallelism = _maxThreadCount;

            var mostAbuIsotopeInternalIndex = _isotopeList.SortedIndexByIntensity[0];
            var mostAbuIsotopeIndex = _isotopeList[mostAbuIsotopeInternalIndex].Index;

            Parallel.ForEach(rows, options, row =>
            {
                Array.Clear(_correlationMap[row], 0, _nScans);
                Array.Clear(_checkedOut[row], 0, _nScans);
                Array.Clear(_distanceMap[row], 0, _nScans);
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

                    var isotopeMzLb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum - 1), charge, isotopeIndex);
                    var isotopeMzUb = (k == 0) ? Ion.GetIsotopeMz(Comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex) : Ion.GetIsotopeMz(Comparer.GetMzAverage(queryMassBinNum + 1), charge, isotopeIndex);
                    
                    if (isotopeMzLb < _run.MinMs1Mz || isotopeMzUb > _run.MaxMs1Mz) continue;
                    var st = _ms1PeakList.BinarySearch(new Ms1Peak(isotopeMzLb, 0));

                    if (st < 0) st = ~st;

                    for (var j = st; j < _ms1PeakList.Count; j++)
                    {
                        var ms1Peak         = _ms1PeakList[j];
                        if (ms1Peak.Mz > isotopeMzUb) break;
                        var col = ms1Peak.Ms1SpecIndex;

                        if (k < 1) // if (k < 4)
                        {
                            if (_featureMatrix[row][col][i] == null || ms1Peak.Intensity > _featureMatrix[row][col][i].Intensity)
                            {
                                _featureMatrix[row][col][i] = ms1Peak;
                                _accurateMass[row][col] = Ion.GetMonoIsotopicMass(ms1Peak.Mz, charge, isotopeIndex);
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
                            var mostAbuPeak = _featureMatrix[row][col][mostAbuIsotopeInternalIndex];
                            if (mostAbuPeak == null) continue;
                            var expectedPeakMz = mostAbuPeak.Mz + (Constants.C13MinusC12*(isotopeIndex - mostAbuIsotopeIndex))/charge;

                            //var mz = Ion.GetIsotopeMz(_accurateMass[row][col], charge, isotopeIndex);
                            if (Math.Abs(expectedPeakMz - ms1Peak.Mz) > MzTolerance.GetToleranceAsTh(ms1Peak.Mz)) continue;

                            // in case of existing isotope peaks, select peaks maximizing envelope similairty
                            if (_featureMatrix[row][col][i] != null)
                            {
                                var tmpPeak = _featureMatrix[row][col][i];
                                var bc1 = _featureMatrix[row][col].GetBhattacharyyaDistance(_isotopeList.EnvelopePdf);
                                _featureMatrix[row][col][i] = ms1Peak;
                                var bc2 = _featureMatrix[row][col].GetBhattacharyyaDistance(_isotopeList.EnvelopePdf);
                                if (bc1 < bc2) _featureMatrix[row][col][i] = tmpPeak;
                            }
                            else
                            {
                                _featureMatrix[row][col][i] = ms1Peak;
                            }
                        }
                    }
                }

                for (var col = 0; col < _nScans; col++)
                {
                    if (!(_accurateMass[row][col] > 0)) continue;
                    var mostAbuPeakIntensity = _featureMatrix[row][col][mostAbuIsotopeInternalIndex].Intensity;
                    var signalToNoiseRatio = mostAbuPeakIntensity / _spectrums[col].MedianIntensity;
                    if (signalToNoiseRatio > 1.4826 && _featureMatrix[row][col].Count(p => p != null && p.Active) >= 3)
                    {
                        _correlationMap[row][col] = _featureMatrix[row][col].GetPearsonCorrelation(_isotopeList.Envelope);
                        _distanceMap[row][col] = _featureMatrix[row][col].GetBhattacharyyaDistance(_isotopeList.EnvelopePdf);
                        if (!observedCharges[row]) observedCharges[row] = true;
                    }
                    else
                    {
                        _accurateMass[row][col] = 0d;
                        //_distanceMap[row][col] = 1.0d;
                    }
                }
            }// end or row for-loop
            );

            var temp = new List<int>();
            for (var i = 0; i < observedCharges.Length; i++) if (observedCharges[i]) temp.Add(i);
            _chargeIndexes = temp.ToArray();
        }
      
        private IEnumerable<ObservedEnvelope> GetSeedCells()
        {
            var seedList = new List<KeyValuePair<double, ObservedEnvelope>>();
            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (!(_accurateMass[i][j] > 0)) continue;

                    var mostAbuPeakIntensity = _featureMatrix[i][j][_isotopeList.SortedIndexByIntensity[0]].Intensity;
                    var signalToNoiseRatio = mostAbuPeakIntensity / _spectrums[j].MedianIntensity;
                    if (signalToNoiseRatio < 3) continue;

                    if (_featureMatrix[i][j].Count(p => p != null && p.Active) < 5) continue;

                    var bcDist = _distanceMap[i][j];
                    if (bcDist > 0.3 && _correlationMap[i][j] < 0.5) continue;
                    var seedCell = new ObservedEnvelope(i, j, _featureMatrix[i][j], _isotopeList);
                    seedList.Add(new KeyValuePair<double, ObservedEnvelope>(bcDist, seedCell));
                }
            }
            return seedList.OrderBy(x => x.Key).Select(x => x.Value);
            //return seedList.OrderByDescending(x => x.Key).Select(x => x.Value);
        }

        private double GetBcDistTh(double normalizedElutionLen)
        {
            if (_queryMass > 15000)
            {
                if (normalizedElutionLen < 0.005) return 1.0;
                if (normalizedElutionLen < 0.01) return 0.6;
                if (normalizedElutionLen < 0.02) return 0.3;
                return 0.2;
            }

            if (normalizedElutionLen < 0.005) return 0.6;
            if (normalizedElutionLen < 0.01) return 0.4;
            if (normalizedElutionLen < 0.02) return 0.2;
            return 0.1;
        }
        private double GetCorrTh(double normalizedElutionLen)
        {
            if (_queryMass > 15000)
            {
                if (normalizedElutionLen < 0.005) return -1;
                if (normalizedElutionLen < 0.01) return 0.001;
                if (normalizedElutionLen < 0.02) return 0.4;
                return 0.6;
            }

            if (normalizedElutionLen < 0.005) return 0.3;
            if (normalizedElutionLen < 0.01) return 0.4;
            if (normalizedElutionLen < 0.02) return 0.6;
            return 0.7;
        }

        private List<Ms1FeatureCluster> FindClusters(bool filtering = false)
        {
            BuildFeatureMatrix(); // should be called first

            var clusters = new List<Ms1FeatureCluster>();
            var tempEnvelope = new double[_isotopeList.Count];
            var tempEnvelope2 = new double[_isotopeList.Count];
            var totalElutionPeriod = _run.GetElutionTime(_ms1ScanNums.Last());

            foreach (var seedCell in GetSeedCells())
            {
                if (_checkedOut[seedCell.Row][seedCell.Col]) continue;
                var seedMass = _accurateMass[seedCell.Row][seedCell.Col];
                var massTol = MzTolerance.GetToleranceAsTh(seedMass);

                var newCluster = new Ms1FeatureCluster(_chargeRange.Min, _ms1ScanNums, _isotopeList);
                Array.Clear(tempEnvelope2, 0, tempEnvelope2.Length);
                seedCell.Peaks.SumEnvelopeTo(tempEnvelope2);
                newCluster.AddMember(seedCell);
                var normalizedElutionLength = (_run.GetElutionTime(newCluster.MaxScanNum) - _run.GetElutionTime(newCluster.MinScanNum)) / totalElutionPeriod;

                var neighbors = new Queue<ObservedEnvelope>();

                neighbors.Enqueue(seedCell); // pick a seed
                _checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    var charge = cell.Row + _chargeRange.Min;

                    var chargeNeighborGap = 1;
                    if (charge > 40) chargeNeighborGap = 12;
                    else if (charge > 30) chargeNeighborGap = 8;
                    else if (charge > 20) chargeNeighborGap = 4;
                    else if (charge > 10) chargeNeighborGap = 2;

                    var minRw = Math.Max(cell.Row - chargeNeighborGap, _chargeIndexes.First());
                    var maxRw = Math.Min(cell.Row + chargeNeighborGap, _chargeIndexes.Last());
                 
                    for (var k = 0; k < 5; k++)
                    {
                        var j = cell.Col;
                        if (k < 3) j += k;
                        else j -= (k - 2);
                        
                        if (j < 0 || j >= _nScans) continue;

                        for (var i = minRw; i <= maxRw; i++)
                        {
                            if (_checkedOut[i][j]) continue;
                            if (!(_accurateMass[i][j] > 0)) continue;
                            if (Math.Abs(seedMass - _accurateMass[i][j]) > massTol) continue;

                            if (_distanceMap[i][j] > GetBcDistTh(normalizedElutionLength) || _correlationMap[i][j] < GetCorrTh(normalizedElutionLength)) continue;
                            Array.Clear(tempEnvelope, 0, _isotopeList.Count);
                            
                            cell.Peaks.SumEnvelopeTo(tempEnvelope);
                            _featureMatrix[i][j].SumEnvelopeTo(tempEnvelope);
                            
                            var newDivergence = _isotopeList.GetBhattacharyyaDistance(tempEnvelope);
                            if (newDivergence < _distanceMap[i][j] && newDivergence < _distanceMap[cell.Row][cell.Col])
                            {
                                var newEnvelope = new ObservedEnvelope(i, j, _featureMatrix[i][j], _isotopeList);
                                neighbors.Enqueue(newEnvelope);
                                newCluster.AddMember(newEnvelope);
                                _checkedOut[i][j] = true;
                                _featureMatrix[i][j].SumEnvelopeTo(tempEnvelope2);
                                normalizedElutionLength = (_run.GetElutionTime(newCluster.MaxScanNum) - _run.GetElutionTime(newCluster.MinScanNum)) / totalElutionPeriod;
                            }
                        }
                    }
                }

                for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) _checkedOut[i][j] = true;

                var summedCorr  = _isotopeList.GetPearsonCorrelation(tempEnvelope2);
                var summedBc = _isotopeList.GetBhattacharyyaDistance(tempEnvelope2);
                if (summedCorr < 0.2 || summedBc > 0.15) continue;
                
                //Console.WriteLine("------------------------{0}", seedMass);
                EvaluateFeature(ref newCluster);

                if (newCluster.Envelopes.Count < 1) continue;
                if (newCluster.GoodEnvelopeCount < 1) continue;
                if (filtering && !newCluster.GoodEnough) continue;

                clusters.Add(newCluster);
            }

            return clusters;
        }

        private void EvaluateFeature(ref Ms1FeatureCluster cluster)
        {
            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;
            var seedMass = _accurateMass[cluster.Envelopes[0].Row][cluster.Envelopes[0].Col];
            var massTol = MzTolerance.GetToleranceAsTh(seedMass);

            cluster.ClearMember();
            for (var col = minCol; col <= maxCol; col++)
            {
                for (var row = minRow; row <= maxRow; row++)
                {
                    var mass = _accurateMass[row][col];
                    if (mass > 0 && Math.Abs(seedMass - mass) < massTol)
                    {
                        var obsEnv = new ObservedEnvelope(row, col, _featureMatrix[row][col], _isotopeList);
                        cluster.AddMember(obsEnv);
                    }
                    else //if (_queryMass > 15000)
                    {
                        var observedPeaks = _spectrums[col].GetAllIsotopePeaks(seedMass, row + _chargeRange.Min, _isotopeList, MzTolerance);
                        var obsEnv = new ObservedEnvelope(row, col, observedPeaks, _isotopeList);
                        if (obsEnv.NumberOfPeaks < 3) continue;
                        cluster.AddMember(obsEnv);    
                    }
                }
            }

            if (cluster.Envelopes.Count < 1) return;
            cluster.NormalizedElutionLength = (_run.GetElutionTime(cluster.MaxScanNum) -
                                               _run.GetElutionTime(cluster.MinScanNum)) /
                                              _run.GetElutionTime(_ms1ScanNums.Last());
            cluster.UpdateScores(_spectrums, seedMass, true);
            if (cluster.GoodEnvelopeCount > 0)
            {
                CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
            }
        }

        private const int MinXicWindowLength = 10;
        
        private void CalculateXicCorrelationOverTimeBetweenIsotopes(Ms1FeatureCluster cluster)
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

            //var n = Math.Min(_isotopeList.Count, 4);
            var xicProfile = new double[_isotopeList.NhighAbundantIsotopes][];
            for (var i = 0; i < xicProfile.Length; i++) xicProfile[i] = new double[colLen];

            var monoMass = cluster.RepresentativeMass;
            var massTol = MzTolerance.GetToleranceAsTh(monoMass);
            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    var mass = _accurateMass[row][col];
                    if (!(mass > 0)) continue;
                    if (Math.Abs(monoMass - mass) > massTol) continue;
                    for (var k = 0; k < xicProfile.Length; k++)
                    {
                        var isoIndex = _isotopeList.SortedIndexByIntensity[k];
                        if (_featureMatrix[row][col][isoIndex] == null) continue;
                        xicProfile[k][col - minCol] += _featureMatrix[row][col][isoIndex].Intensity;
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

        /*
        public double GetChiSquareSignificanceTest(double[] observedIsotopeEnvelop)
        {
            var k = Count - 1;
            var x = 0d;
            var s2 = 0d;

            for (var i = 0; i < Count; i++)
            {
                s2 += observedIsotopeEnvelop[i];
            }
            for (var i = 0; i < Count; i++)
            {
                var p = EnvelopePdf[i];
                var q = observedIsotopeEnvelop[i]/s2;
                x += ((p - q) * (p - q)) / (p + q);
            }

            x *= 0.5;
            var pvalue = ChiSquared.CDF(k, x);

            return pvalue;
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
        */
    }

}