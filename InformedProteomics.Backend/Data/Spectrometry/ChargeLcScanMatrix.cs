using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {
        public static void Initilize(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int numBits = 27)
        {
            _run = run;
            _minScanCharge = minScanCharge;
            _maxScanCharge = maxScanCharge;

            _mzTolerance = new Tolerance(10);
            _comparer = new MzComparerWithBinning(numBits);

            _ms1ScanNums = run.GetMs1ScanVector();
            _ms1ScanNumToIndexMap = new int[_ms1ScanNums.Last() - _ms1ScanNums.First() + 1];

            _nScans = _ms1ScanNums.Length;
            
            for (var i = 0; i < _ms1ScanNums.Length; i++)
                _ms1ScanNumToIndexMap[_ms1ScanNums[i] - _ms1ScanNums[0]] = i;

            _ms1PeakList = new List<LcMsPeak>();
            _cachedSpectrum = new List<ChargeLcScanSpectrum>();

            foreach (var scanNum in _ms1ScanNums)
            {
                var spec = new ChargeLcScanSpectrum(_run.GetSpectrum(scanNum));
                _cachedSpectrum.Add(spec);
                foreach (var peak in spec.Peaks)
                    _ms1PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum));
            }
            _ms1PeakList.Sort();
        }

        public ChargeLcScanMatrix(LcMsRun run = null, int minScanCharge = 2, int maxScanCharge = 60, int numBits = 27)
        {
            if (run != null) Initilize(run, minScanCharge, maxScanCharge, numBits);
            if (_run == null) throw new Exception("ChargeLcScanMatrix has not been initialized");

            _intensityMapFull = new double[MaxChargeLength][];
            _mostAbuIsotopeMzMap = new double[MaxChargeLength][];
            _correlationMap = new double[MaxChargeLength][];
            _featureMatrix = new double[MaxChargeLength][][];
            _checkedOut = new bool[MaxChargeLength][];
            _isotopeMzLb = new double[MaxChargeLength][];
            _isotopeMzUb = new double[MaxChargeLength][];

            for (var i = 0; i < MaxChargeLength; i++)
            {
                _intensityMapFull[i] = new double[_nScans];
                _mostAbuIsotopeMzMap[i] = new double[_nScans];
                _checkedOut[i] = new bool[_nScans];
                _correlationMap[i] = new double[_nScans];
                _featureMatrix[i] = new double[_nScans][];
                _isotopeMzLb[i] = new double[MaxEnvelopeLength];
                _isotopeMzUb[i] = new double[MaxEnvelopeLength];

                for (var j = 0; j < _nScans; j++)
                    _featureMatrix[i][j] = new double[MaxEnvelopeLength];
            }
        }
        
        public static ChargeLcScanMatrix GetInstance()
        {
            return new ChargeLcScanMatrix();
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
        private static int[] _ms1ScanNumToIndexMap;
        private static MzComparerWithBinning _comparer;
        private static Tolerance _mzTolerance;
        private static List<ChargeLcScanSpectrum> _cachedSpectrum;
        private static List<LcMsPeak> _ms1PeakList;
        private const int MaxEnvelopeLength = 30;
        private const int MaxChargeLength = 40;
        private const double CorrLowerBound = 0.5d;
        //private static SavitzkyGolaySmoother _smoother = new SavitzkyGolaySmoother(9, 2);

        private readonly double[][] _intensityMapFull;
        private readonly double[][] _mostAbuIsotopeMzMap;
        private readonly double[][] _correlationMap;
        private readonly double[][][] _featureMatrix;
        private readonly bool[][] _checkedOut;

        private readonly double[][] _isotopeMzLb;
        private readonly double[][] _isotopeMzUb;

        private IsotopeList _isotopeList;
        private int[] _chargeIndexes;
        private IntRange _chargeRange;
        private int _minChargeRangeLength;
        private double _queryMass;
        

        private void SetQueryMass(double queryMass)
        {
            _queryMass = queryMass;
            _chargeRange = GetScanChargeRange(_queryMass);

            _minChargeRangeLength = Math.Min((int) Math.Ceiling((queryMass/1000d) - 15), 30);
            _minChargeRangeLength = Math.Max(_minChargeRangeLength, 0);
        }

        private IntRange GetScanChargeRange(double mass)
        {
            if (mass < 5000.0d) return new IntRange(_minScanCharge, 10);
            
            var chargeLb = (int)Math.Max(_minScanCharge, Math.Floor((13.0/2.5) * (mass / 10000d) - 0.6));
            var chargeUb = (int)Math.Min(_maxScanCharge, Math.Ceiling(16 * (mass / 10000d) + 5));
            if (chargeUb - chargeLb + 1 > MaxChargeLength) chargeUb = chargeLb + MaxChargeLength - 1;

            return new IntRange(chargeLb, chargeUb);
        }

        private void BuildFeatureMatrix()
        {
            var queryMassBinNum = _comparer.GetBinNumber(_queryMass);
            _isotopeList = new IsotopeList(_comparer.GetMzAverage(queryMassBinNum), MaxEnvelopeLength);

            for (var row = 0; row < _chargeRange.Length; row++)
            {
                Array.Clear(_intensityMapFull[row], 0, _nScans);
                Array.Clear(_mostAbuIsotopeMzMap[row], 0, _nScans);
                Array.Clear(_correlationMap[row], 0, _nScans);
                Array.Clear(_checkedOut[row], 0, _nScans);

                for (var col = 0; col < _nScans; col++)
                {
                    Array.Clear(_featureMatrix[row][col], 0, _featureMatrix[row][col].Length);
                }
            }
            
            for (var row = 0; row < _chargeRange.Length; row++)
            {
                var charge = row + _chargeRange.Min;
                for (var i = 0; i < _isotopeList.Count; i++)
                {
                    var isotopeIndex = _isotopeList[i].Index;
                    var mzMin = Ion.GetIsotopeMz(_comparer.GetMzStart(queryMassBinNum), charge, isotopeIndex);
                    var mzMax = Ion.GetIsotopeMz(_comparer.GetMzEnd(queryMassBinNum), charge, isotopeIndex);
                    _isotopeMzLb[row][i] = mzMin - _mzTolerance.GetToleranceAsTh(mzMin);
                    _isotopeMzUb[row][i] = mzMax + _mzTolerance.GetToleranceAsTh(mzMax);
                }
            }

            var highestPeakIndex        = new int[_nScans];
            var highestPeakIntensity    = new double[_nScans];
            var observedCharges         = new bool[_chargeRange.Length];
            // Build FeatureMatrix...Starting from the most abundant isotope
            for (var row = 0; row < _chargeRange.Length; row++)
            {
                Array.Clear(highestPeakIndex, 0, highestPeakIndex.Length);
                Array.Clear(highestPeakIntensity, 0, highestPeakIntensity.Length);

                for (var k = 0; k < _isotopeList.Count; k++)
                {
                    var i = _isotopeList.SortedIndexByIntensity[k]; // internal isotope index

                    if (_isotopeMzLb[row][i] < _run.MinMs1Mz || _isotopeMzUb[row][i] > _run.MaxMs1Mz) continue;

                    var st = _ms1PeakList.BinarySearch(new LcMsPeak(_isotopeMzLb[row][i], 0, 0));
                    if (st < 0) st = ~st;

                    for (var j = st; j < _ms1PeakList.Count; j++)
                    {
                        var ms1Peak = _ms1PeakList[j];
                        if (ms1Peak.Mz > _isotopeMzUb[row][i]) break;
                        var col = _ms1ScanNumToIndexMap[ms1Peak.ScanNum - _ms1ScanNums.First()];

                        if (k <= 3)
                        {
                            _featureMatrix[row][col][i] = Math.Max(ms1Peak.Intensity, _featureMatrix[row][col][i]);
                            if (_featureMatrix[row][col][i] > highestPeakIntensity[col])
                            {
                                highestPeakIndex[col] = i;
                                highestPeakIntensity[col] = _featureMatrix[row][col][i];

                                if (k == 0) _mostAbuIsotopeMzMap[row][col] = ms1Peak.Mz;
                            }
                        }
                        else
                        {
                            if (_featureMatrix[row][col][i] > 0)
                            {
                                // if existing peak matches better, then skip
                                if (highestPeakIntensity[col] > 0)
                                {
                                    var expectedIntensity = (highestPeakIntensity[col] * _isotopeList[i].Ratio) / _isotopeList[highestPeakIndex[col]].Ratio;
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
            }

            var temp = new List<int>();
            for (var i = 0; i < observedCharges.Length; i++) if (observedCharges[i]) temp.Add(i);
            _chargeIndexes = temp.ToArray();

            // Fill _correlationMap array
            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    _intensityMapFull[i][j] = _featureMatrix[i][j].Sum();
                    if (!(_intensityMapFull[i][j] > 0)) continue;
                    _correlationMap[i][j] = _isotopeList.GetPearsonCorrelation(_featureMatrix[i][j]);
                }
            }
        }


        private List<ChargeLcScanCluster> FindClusters(double clusteringScoreCutoff = 0.7)
        {
            BuildFeatureMatrix(); // should be called first

            var clusters = new List<ChargeLcScanCluster>();
            var tempEnvelope = new double[_isotopeList.Count];

            foreach (var seedCell in GetSeedCells(CorrLowerBound))
            {
                if (_checkedOut[seedCell.Row][seedCell.Col]) continue;

                var chargeNeighborGap = 1;
                var scanNeighborGap = 1;
                if (seedCell.Row + _chargeRange.Min >= 20)
                {
                    scanNeighborGap = 2;
                    chargeNeighborGap = 2;
                }

                var seedScore = _correlationMap[seedCell.Row][seedCell.Col];
                var newCluster = new ChargeLcScanCluster(_chargeRange.Min, _ms1ScanNums, _isotopeList);
                newCluster.AddMember(seedCell, _featureMatrix[seedCell.Row][seedCell.Col], _correlationMap[seedCell.Row][seedCell.Col]);

                var neighbors = new Queue<ChargeLcScanCell>();

                neighbors.Enqueue(seedCell); // pick a seed
                _checkedOut[seedCell.Row][seedCell.Col] = true;
                
                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();

                    for (var k = Math.Max(cell.Row - chargeNeighborGap, _chargeIndexes.First()); k <= Math.Min(cell.Row + chargeNeighborGap, _chargeIndexes.Last()); k++)
                    {
                        for (var l = Math.Max(cell.Col - scanNeighborGap, 0); l <= Math.Min(cell.Col + scanNeighborGap, _nScans - 1); l++)
                        {
                            var distFromSeed = Math.Abs(seedCell.Col - l);

                            if (_checkedOut[k][l] || (_correlationMap[k][l] < CorrLowerBound && distFromSeed > 1) ||
                                (_correlationMap[k][l] < 0.2 && distFromSeed <= 1)) continue;

                            for (var t = 0; t < tempEnvelope.Length; t++) tempEnvelope[t] = newCluster.ClusteringEnvelope[t] + _featureMatrix[k][l][t];
                            var newScore = _isotopeList.GetPearsonCorrelation(tempEnvelope);

                            var cellCorrTh = Math.Min(0.6 + (distFromSeed - 1) * 0.03, 0.9);
                            var summmedCorrTh = Math.Min(0.8 + (distFromSeed - 1) * 0.015, 0.95);

                            if (!(newCluster.ClusteringScore < newScore ||
                                (distFromSeed < 5 && seedScore < newScore) ||
                                (_correlationMap[k][l] > cellCorrTh && newScore > summmedCorrTh))) continue;

                            var newMember = new ChargeLcScanCell(k, l);
                            neighbors.Enqueue(newMember);
                            newCluster.AddMember(newMember, _featureMatrix[k][l], newScore);
                            _checkedOut[k][l] = true;
                        }
                    }
                }

                for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) _checkedOut[i][j] = true;

                // high charges should be spread over a range
                
                if (newCluster.ClusteringScore < clusteringScoreCutoff || newCluster.ChargeLength < _minChargeRangeLength) continue;
                if (newCluster.ScanLength <= 1) continue;

                clusters.Add(newCluster);
            }

            return clusters;
        }


        private List<ChargeLcScanCluster> GetFilteredClusters(IEnumerable<ChargeLcScanCluster> clusters, bool filtering = true)
        {
            var filteredClusters = new List<ChargeLcScanCluster>();

            foreach (var cluster in clusters)
            {
                GetAccurateMassInfo(cluster);

                if (cluster.Members.Count <= 1 || cluster.ScanLength <= 1) continue;
                if (filtering && cluster.ChargeLength < _minChargeRangeLength) continue;
                if (filtering && cluster.GetScore(ChargeLcScanScore.RankSum) + cluster.GetScore(ChargeLcScanScore.RankSumSummed) < 10) continue;
                if (filtering && cluster.GetScore(ChargeLcScanScore.HyperGeometric) + cluster.GetScore(ChargeLcScanScore.HyperGeometricSummed) < 60) continue;
                if (filtering && cluster.GetScore(ChargeLcScanScore.EnvelopeCorrelation) < 0.5) continue;
                if (filtering && cluster.GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.6) continue;

                filteredClusters.Add(cluster);
            }

            return filteredClusters;
        }
        
        private void GetAccurateMassInfo(ChargeLcScanCluster cluster)
        {
            var bestEnvelopeCorrelation = -1.0d;
            var bestRankSumScore = 1.0d;
            var bestHyperGeometricScore = 1.0d;
            
            var summedK = 0;
            var summedK1 = 0;
            var summedN = 0;
            var summedN1 = 0;

            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;

            var memberCorr = new List<double>();
            var memberMass = new List<double>();
            var memberRankSum = new List<double>();
            
            cluster.ClearMember();

            for (var col = minCol; col <= maxCol; col++)
            {
                var scanNum = _ms1ScanNums[col];
                for (var row = minRow; row <= maxRow; row++)
                {
                    var charge = row + _chargeRange.Min;
                    var mostAbuIsoMz = _mostAbuIsotopeMzMap[row][col];
                    
                    if (!(mostAbuIsoMz > 0)) continue;

                    var envelope = new double[_isotopeList.Count];
                    var mass = Ion.GetMonoIsotopicMass(mostAbuIsoMz, charge, _isotopeList.GetMostAbundantIsotope().Index);
                    var spectrum = GetSpectrum(col);
                    Peak[] peaks;
                    Peak[] observedPeaks;
                    int[] observedPeakRanks;

                    spectrum.GetRankSum(_isotopeList, charge, mass, _mzTolerance, out peaks, out observedPeaks, out observedPeakRanks);
                    
                    var k1 = 0;
                    for (var i = 0; i < _isotopeList.Count; i++)
                    {
                        if (observedPeaks[i] == null) continue;
                        envelope[i] = observedPeaks[i].Intensity;
                        k1++;
                    }
                    
                    var envCorr = _isotopeList.GetPearsonCorrelation(envelope);

                    var binNum = spectrum.RankingWindowTolerance.GetBinNumber(mostAbuIsoMz);
                    var maxMz = spectrum.RankingWindowTolerance.GetMzEnd(binNum);
                    var minMz = spectrum.RankingWindowTolerance.GetMzStart(binNum);
                    var n = (int)Math.Round((maxMz - minMz) / (0.5*(_mzTolerance.GetToleranceAsTh(maxMz)+_mzTolerance.GetToleranceAsTh(minMz))));
                    var k = peaks.Length;
                    var n1 = _isotopeList.Count;
                    
                    if (envCorr < CorrLowerBound || k1 < 3) continue;

                    var ranksumScore = CalculateRankSumScore(peaks, observedPeaks, observedPeakRanks);
                    var hyperGeoScore = CalculateHyperGeometricScore(n, k, n1, k1);

                    //if (ranksumScore < 8 || hyperGeoScore < 10) continue;
                    
                    if (bestHyperGeometricScore < hyperGeoScore) bestHyperGeometricScore = hyperGeoScore;
                    if (bestRankSumScore < ranksumScore) bestRankSumScore = ranksumScore;
                    if (bestEnvelopeCorrelation < envCorr) bestEnvelopeCorrelation = envCorr;

                    summedN += n;
                    summedK += k;
                    summedK1 += k1;
                    summedN1 = n1;
                    memberRankSum.Add(ranksumScore);

                    cluster.AddMember(new ChargeLcScanCell(row, col), envelope);
                    memberCorr.Add(envCorr);
                    memberMass.Add(mass);
                }
            }
            
            cluster.SetScore(ChargeLcScanScore.EnvelopeCorrelation, bestEnvelopeCorrelation);
            cluster.SetScore(ChargeLcScanScore.RankSum, bestRankSumScore);
            cluster.SetScore(ChargeLcScanScore.HyperGeometric, bestHyperGeometricScore);

            var summedRankSum = 0d;
            var summedHyperGeo = 0d;
            var summedCorr = 0d;
            if (cluster.Members.Count > 0)
            {
                memberRankSum.Sort();
                memberRankSum.Reverse();
                for (var k = 0; k < memberRankSum.Count && k < 10; k++)
                {
                    summedRankSum += memberRankSum[k];
                }
                summedRankSum /= Math.Min(memberRankSum.Count, 10);
                summedHyperGeo = CalculateHyperGeometricScore(summedN, summedK, summedN1, summedK1);
                summedCorr = GetBestSummedCorrelation(cluster, memberCorr, memberMass);
            }

            cluster.SetScore(ChargeLcScanScore.EnvelopeCorrelationSummed, summedCorr);
            cluster.SetScore(ChargeLcScanScore.RankSumSummed, summedRankSum);
            cluster.SetScore(ChargeLcScanScore.HyperGeometricSummed, summedHyperGeo);
        }

        private double GetBestSummedCorrelation(ChargeLcScanCluster cluster, List<double> cellCorr, List<double> cellMass)
        {
            var cellCorrTemp = new double[cellCorr.Count];
            var index = Enumerable.Range(0, cellCorrTemp.Length).ToArray();

            Array.Copy(cellCorr.ToArray(), cellCorrTemp, cellCorr.Count);
            Array.Sort(cellCorrTemp, index);
            Array.Reverse(index);

            var summedEnvelop = new double[_isotopeList.Count];
            var summedCorr = 0d;
            var tempMass = new List<KeyValuePair<int, double>>();

            foreach (var j in index)
            {
                var tempEnvelop = new double[_isotopeList.Count];
                for (var k = 0; k < _isotopeList.Count; k++) tempEnvelop[k] = summedEnvelop[k] + cluster.MemberEnvelope[j][k];
                var tempCorr = _isotopeList.GetPearsonCorrelation(tempEnvelop);

                if (!(tempCorr > summedCorr)) continue;

                summedEnvelop = tempEnvelop;
                summedCorr = tempCorr;

                if (tempMass.Count < 10) tempMass.Add(new KeyValuePair<int, double>(j, cellMass[j]));
            }
            cluster.SummedEnvelope = summedEnvelop;

            var medIdx = (int) Math.Round(tempMass.Count*0.5);
            var i = 0;
            foreach (var m in tempMass.OrderBy(x => x.Value))
            {
                if (i == medIdx)
                {
                    var row = cluster.Members[m.Key].Row;
                    var col = cluster.Members[m.Key].Col;
                    cluster.SetRepresentativeMass(m.Value, _mostAbuIsotopeMzMap[row][col], row + _chargeRange.Min, _ms1ScanNums[col]);
                    break;
                }
                i++;
            }
            return summedCorr;
        }

        private static double CalculateHyperGeometricScore(int n, int k, int n1, int k1)
        {
            var pvalue = FitScoreCalculator.GetHyperGeometricPvalue(n, k, n1, k1);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        private double CalculateRankSumScore(Peak[] peakList, Peak[] isotopePeaks, int[] isotopePeakRanks)
        {
            var n = 0;
            var r = 0d;

            for(var k = 0; k < _isotopeList.Count; k++)
            {
                var internalIndex = _isotopeList.SortedIndexByIntensity[k];

                if (_isotopeList[internalIndex].Ratio < 0.4) break;

                if (isotopePeakRanks[internalIndex] > 0)
                {
                    n++;
                    r += isotopePeakRanks[internalIndex];
                }
            }
            
            if (peakList.Length == n) return 0d;
            var pvalue = FitScoreCalculator.GetRankSumPvalue(peakList.Length, n, r);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }
      
       

        private IEnumerable<ChargeLcScanCell> GetSeedCells(double seedCutoff)
        {
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();
            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (_correlationMap[i][j] < seedCutoff) continue;
                    seedCells.Add(new KeyValuePair<double, ChargeLcScanCell>(_correlationMap[i][j], new ChargeLcScanCell(i, j)));
                }
            }

            return seedCells.OrderByDescending(x => x.Key).Select(x => x.Value);
        }



        private ChargeLcScanSpectrum GetSpectrum(int col)
        {
            return _cachedSpectrum[col];
        }

        public double[][] GetCorrelationMap()
        {
            return _correlationMap;
        }
        
        public double[][] GetIntensityMap()
        {
            return _intensityMapFull;
        }

        
    }

    public class IsotopeList : List<Isotope>
    {
        public double MonoIsotopeMass { get; private set; }
        private readonly double[] _envelope;
        //private readonly int[] _sortedIndex;
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

            for (var i = 0; i < Count; i++)
            {
                var rankingIndex = isotopeRankings[this[i].Index] - 1;
                SortedIndexByIntensity[rankingIndex] = i;
            }
        }

        public Isotope GetIsotopeRankedAt(int ranking)
        {
            return this[SortedIndexByIntensity[ranking - 1]];
        }

        public double GetPearsonCorrelation(double[] observedIsotopeEnvelop)
        {
            return FitScoreCalculator.GetPearsonCorrelation(_envelope, observedIsotopeEnvelop, _envelope.Length);
        }

        public Isotope GetMostAbundantIsotope()
        {
            return GetIsotopeRankedAt(1);
        }
    }


    class ChargeLcScanSpectrum : Spectrum
    {
        internal ChargeLcScanSpectrum(Spectrum spec)
            : base(spec.Peaks, spec.ScanNum)
        {
            _comparer = new MzComparerWithBinning(17);
            _minBinNum = _comparer.GetBinNumber(MinMz);
            _maxBinNum = _comparer.GetBinNumber(MaxMz);
            _peakIndex = new IntRange[3][];
            _peakRanking = new int[3][][];
            for (var i = 0; i < 3; i++)
            {
                _peakIndex[i] = new IntRange[_maxBinNum - _minBinNum + 1];
                _peakRanking[i] = new int[_maxBinNum - _minBinNum + 1][];
            }


            for (var binNum = _minBinNum; binNum <= _maxBinNum; binNum++)
            {
                var minMz = _comparer.GetMzAverage(binNum - 1);
                var maxMz = _comparer.GetMzAverage(binNum + 1);

                var startIndex = Array.BinarySearch(Peaks, new Peak(minMz - float.Epsilon, 0));
                if (startIndex < 0) startIndex = ~startIndex;

                // go up
                var intensities = new List<double>();
                var intensities1 = new List<double>();
                var intensities2 = new List<double>();
                var minMz0 = _comparer.GetMzStart(binNum);
                var maxMz0 = _comparer.GetMzEnd(binNum);

                var minMz1 = _comparer.GetMzAverage(binNum - 1);
                var maxMz1 = _comparer.GetMzAverage(binNum);
                var minMz2 = maxMz1;
                var maxMz2 = _comparer.GetMzAverage(binNum + 1);

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
                        intensities.Add(Peaks[i].Intensity);
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
                    _peakRanking[0][binNum - _minBinNum] = ArrayUtil.GetRankings(intensities);
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

        internal MzComparerWithBinning RankingWindowTolerance { get { return _comparer; } }

        internal int GetRanks(int binShift, int binIndex, int globalIndex)
        {
            int ret = 0;

            if (globalIndex >= _peakIndex[binShift][binIndex].Min && globalIndex <= _peakIndex[binShift][binIndex].Max)
                ret = _peakRanking[binShift][binIndex][globalIndex - _peakIndex[binShift][binIndex].Min];

            return ret;
        }

        internal void GetRankSum(IsotopeList isotopeList, int charge, double monoIsotopeMass, Tolerance tolerance, out Peak[] peaks, out Peak[] observedPeaks, out int[] observedPeakRanks)
        {
            var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeList.GetMostAbundantIsotope().Index);
            var mostAbundantIsotopePeakIndex = FindPeakIndex(mostAbundantIsotopeMz, tolerance);
            if (mostAbundantIsotopePeakIndex < 0) throw new ArgumentOutOfRangeException();

            var binNum = _comparer.GetBinNumber(mostAbundantIsotopeMz);
            var binIndex = binNum - _minBinNum;
            var binShift = 0;

            var d0 = Math.Abs(_comparer.GetMzAverage(binNum) - mostAbundantIsotopeMz);
            var d1 = Math.Abs(_comparer.GetMzStart(binNum) - mostAbundantIsotopeMz);
            var d2 = Math.Abs(_comparer.GetMzEnd(binNum) - mostAbundantIsotopeMz);

            if (d1 < d2 && d1 < d0) binShift = 1;
            else if (d2 < d1 && d2 < d0) binShift = 2;

            if (_peakIndex[binShift][binIndex] == null) throw new ArgumentOutOfRangeException();

            peaks = new Peak[_peakIndex[binShift][binIndex].Max - _peakIndex[binShift][binIndex].Min + 1];
            Array.Copy(Peaks, _peakIndex[binShift][binIndex].Min, peaks, 0, peaks.Length);

            var mostAbuInternalIndex = isotopeList.SortedIndexByIntensity[0];
            observedPeaks = new Peak[isotopeList.Count];
            observedPeakRanks = new int[isotopeList.Count];
            observedPeaks[mostAbuInternalIndex] = Peaks[mostAbundantIsotopePeakIndex];
            observedPeakRanks[mostAbuInternalIndex] = GetRanks(binShift, binIndex, mostAbundantIsotopePeakIndex);

            // go down
            var peakIndex = mostAbundantIsotopePeakIndex - 1;
            for (var isotopeInternalIndex = isotopeList.SortedIndexByIntensity[0] - 1; isotopeInternalIndex >= 0; isotopeInternalIndex--)
            {
                var isotopeIndex = isotopeList[isotopeInternalIndex].Index; //selectedIsotopeIndex[isotopeInternalIndex];
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
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
                        }
                    }
                }
            }
        }

        internal double MinMz { get { return Peaks[0].Mz; } }
        internal double MaxMz { get { return Peaks[Peaks.Length - 1].Mz; } }

        private readonly MzComparerWithBinning _comparer;
        private readonly IntRange[][] _peakIndex;
        private readonly int[][][] _peakRanking;
        private readonly int _minBinNum;
        private readonly int _maxBinNum;
    }
}
