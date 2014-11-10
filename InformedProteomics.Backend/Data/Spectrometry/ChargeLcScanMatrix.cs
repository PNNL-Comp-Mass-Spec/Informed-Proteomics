using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing.Design;
using System.Drawing.Text;
using System.Linq;
using System.Net.Sockets;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;
using MultiDimensionalPeakFinding;
using MathNet.Numerics.Distributions;


namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {
        public ChargeLcScanMatrix(LcMsRun run, int minCharge = 2, int maxCharge = 50, int numBits = 27, bool useCache = true)
        {
            _minCharge = minCharge;
            _maxCharge = maxCharge;

            _ms1ScanNums = run.GetMs1ScanVector();
            _ms1ScanNumToIndexMap = new Dictionary<int, int>();
            for (var i = 0; i < _ms1ScanNums.Length; i++)
            {
                _ms1ScanNumToIndexMap.Add(_ms1ScanNums[i], i);
            }
            

            _comparer = new MzComparerWithBinning(numBits);

            _nScans = _ms1ScanNums.Length;
            _nCharges = _maxCharge - _minCharge + 1;

            _run = run;
            _clusterMap = new Dictionary<int, IEnumerable<ChargeLcScanCluster>>();
            _binNumToMs2ScanNumsMap = new Dictionary<int, IEnumerable<int>>();

            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);

            if (useCache)
            {
                //Cache Xic data
                _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];
                for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
                {
                    var mzStart = _comparer.GetMzStart(binNum);
                    var mzEnd = _comparer.GetMzEnd(binNum);
                    var mzHalf = 0.5*(mzEnd - mzStart);
                    var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart - mzHalf, mzEnd + mzHalf);
                    _cachedXic[binNum - _cachedMinBinNum] = xic;
                }

                // Generate Ranking information 
                _cachedRanking = new int[_nScans][];
                _cachedNumPeaks = new int[_nScans];
                var tempXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1];
                
                for (var i = 0; i < _nScans; i++)
                {
                    for (var j = 0; j < _cachedMaxBinNum - _cachedMinBinNum + 1; j++)
                    {
                        tempXic[j] = _cachedXic[j][i];
                    }
                    _cachedNumPeaks[i] = tempXic.Where(x => x > 0).Count();
                    _cachedRanking[i] = ArrayUtil.GetRankings(tempXic, 0.1d);
                }
            }

            _intensityMapFull           = new double[_nCharges][];
            _intensityMapMostAbundant = new double[_nCharges][];
            _correlationMap     = new double[_nCharges][];
            _xicMatrix          = new double[_nCharges][][];
            _rankingMatrix      = new int[_nCharges][][];

            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMapFull[i]            = new double[_nScans];
                _intensityMapMostAbundant[i]    = new double[_nScans];

                _correlationMap[i] = new double[_nScans];
                _xicMatrix[i] = new double[_nScans][];
                _rankingMatrix[i] = new int[_nScans][];

                for (var j = 0; j < _nScans; j++)
                {
                    _xicMatrix[i][j] = new double[MaxEnvelopeLength];
                    _rankingMatrix[i][j] = new int[MaxEnvelopeLength];
                }
            }

            _observedCharges = new bool[_nCharges];
        }

        public IEnumerable<ChargeScanRange> GetProbableChargeScanRegions(double monoIsotopicMass)
        {
            var clusters = GetProbableChargeScanClusters(monoIsotopicMass);

            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            IEnumerable<int> ms2Scans;

            if (_binNumToMs2ScanNumsMap.TryGetValue(binNumber, out ms2Scans)) return ms2Scans;

            var clusters = GetProbableChargeScanClusters(monoIsotopicMass);
            var ms2ScanSet = new HashSet<int>();

            foreach (var cluster in clusters) ms2ScanSet.UnionWith(GetMs2ScanNumbers(cluster));
            
            var ms2ScanList = ms2ScanSet.ToList();
            ms2ScanList.Sort();
            
            return _binNumToMs2ScanNumsMap[binNumber] = ms2ScanList;
        }

        public IEnumerable<ChargeScanRange> GetAllClusters(double monoIsotopicMass, out double[][] scores, out HashSet<int>[] ms2ScanList)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            _monoIsotopicMass = _comparer.GetMzAverage(binNumber);
            //_monoIsotopicMass = monoIsotopicMass;

            //var clusters = FindClusters(3.3, 50, 4);
            var clusters = FindClusters(2.0, 30, 2);

            scores = new double[6][];
            for (var i = 0; i < scores.Length; i++) scores[i] = new double[clusters.Count];

            ms2ScanList = new HashSet<int>[clusters.Count];

            for(var i = 0; i < clusters.Count; i++)
            {
                scores[0][i] = clusters[i].CorrBetweenObsTh;
                scores[1][i] = CalculateXicCorrelationOverTimeBetweenIsotopes(clusters[i]);
                scores[2][i] = CalculateXicCorrelationOverTimeBetweenCharges(clusters[i]);

                scores[3][i] = CalculateHyperGeometricScore(clusters[i]);
                scores[4][i] = CalculateRankSumScore(clusters[i]);

                ms2ScanList[i] = GetMs2ScanNumbers(clusters[i]);
            }
            
            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }

        public MzComparerWithBinning GetMzComparerWithBinning()
        {
            return _comparer;
        }

        private readonly LcMsRun _run;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly int _nCharges;
        private readonly int _nScans;
        private readonly int[] _ms1ScanNums;

        private readonly double[][] _intensityMapFull;
        private readonly double[][] _intensityMapMostAbundant;
        private readonly double[][] _correlationMap;

        private readonly double[][][] _xicMatrix;
        private readonly int[][][] _rankingMatrix;

        private readonly bool[] _observedCharges;

        private double _monoIsotopicMass;
        private double[] _envelope;
        private int _mostAbundantIsotopeIndex;

        private int[] _envelopeIndex = null;
        private const int MaxEnvelopeLength = 50;
        private int[] _chargeIndexes;
        private readonly MzComparerWithBinning _comparer;

        // caching
        private readonly double[][] _cachedXic;
        private readonly int _cachedMinBinNum;
        private readonly int _cachedMaxBinNum;
        private readonly int[][] _cachedRanking;
        private readonly int[] _cachedNumPeaks;

        private readonly Dictionary<int, IEnumerable<ChargeLcScanCluster>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;


        private readonly Dictionary<int, int> _ms1ScanNumToIndexMap;
        
        
        internal IEnumerable<ChargeLcScanCluster> GetProbableChargeScanClusters(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            IEnumerable<ChargeLcScanCluster> clusters;

            if (_clusterMap.TryGetValue(binNumber, out clusters))
            {
                return clusters;
            }
            
            _monoIsotopicMass = _comparer.GetMzAverage(binNumber);

            clusters = FindClusters(3.3, 50, 4);
            return _clusterMap[binNumber] = clusters;
        }
        
        internal IEnumerable<ChargeLcScanCluster> GetMs1Features(double monoIsotopicMass)
        {
            _monoIsotopicMass = monoIsotopicMass;
            var clusters = FindClusters(3.5, 50, 4);

            foreach (var c in clusters)
            {
                c.IntensityAlongCol = new double[c.MaxCol - c.MinCol + 1];
                for (var col = c.MinCol; col <= c.MaxCol; col++)
                {
                    for (var row = c.MinRow; row <= c.MaxRow; row++)
                    {
                        c.IntensityAlongCol[col - c.MinCol] += _intensityMapFull[row][col];
                    }
                }
            }

            return clusters;
        }

        //private const int MaxMsScanNumPerCluster = 10;
        private HashSet<int> GetMs2ScanNumbers(ChargeLcScanCluster cluster)
        {
            var result = new HashSet<int>();
            //var ms2ScanNumbers = new List<KeyValuePair<int, double>>();
            
            for (var row = cluster.MinRow; row <= cluster.MaxRow; row++)
            {
                var charge = _minCharge + row;
                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(_monoIsotopicMass, charge, _envelopeIndex[_mostAbundantIsotopeIndex]);
                var ms2 = _run.GetFragmentationSpectraScanNums(mostAbundantIsotopeMz);

                //var ms2Set1 = ms2.Where(scanNum => scanNum >= _ms1ScanNums[cluster.MinCol] && scanNum <= _ms1ScanNums[cluster.MaxCol]); //.ToArray();
                //for (var col = cluster.MinCol; col < cluster.MaxCol; col++)
                var maxAbundance = 0d;
                var selectedMs2ScanNum = 0;
                foreach(var ms2ScanNum in ms2)
                {
                    if (ms2ScanNum < _ms1ScanNums[cluster.MinCol] || ms2ScanNum > _ms1ScanNums[cluster.MaxCol])
                        continue;
                    
                    var prevMs1ScanNum = _run.GetPrevScanNum(ms2ScanNum, 1);
                    var nextMs1ScanNum = _run.GetNextScanNum(ms2ScanNum, 1);
                    var col1 = _ms1ScanNumToIndexMap[prevMs1ScanNum];
                    var col2 = _ms1ScanNumToIndexMap[nextMs1ScanNum];

                    var abundance = _intensityMapMostAbundant[row][col1] + _intensityMapMostAbundant[row][col2];

                    if (abundance > maxAbundance)
                    {
                        maxAbundance = abundance;
                        selectedMs2ScanNum = ms2ScanNum;
                    }
                    //ms2ScanNumbers.Add(new KeyValuePair<int, double>(ms2ScanNum, abundance));
                }
                //ms2ScanNumbers.UnionWith(ms2Set1);
                /*
                if (ms2Set1.Length > 0) continue;
                
                var minCol = Math.Max(cluster.MinCol - 1, 0);
                var maxCol = Math.Min(cluster.MaxCol + 1, _ms1ScanNums.Length - 1);

                var ms2Set2 = ms2.Where(scanNum => scanNum >= _ms1ScanNums[minCol] && scanNum <= _ms1ScanNums[maxCol]).ToArray();
                ms2ScanNumbers.UnionWith(ms2Set2);*/
                if (selectedMs2ScanNum > 0) result.Add(selectedMs2ScanNum);
            }
            
            /*
            foreach (var ms2IntensityPair in ms2ScanNumbers.OrderByDescending(x => x.Value))
            {
                result.Add(ms2IntensityPair.Key);
                if (result.Count > MaxMsScanNumPerCluster) break;
            }*/

            return result;
        }
        
        private double CalculateRankSumScore(ChargeLcScanCluster cluster)
        {
            if (_cachedXic == null)
            {
                return 0;
            }

            var minPvalue = 1.0d;

            foreach (var cell in cluster.Members)
            {
                var r1 = 0;
                var n1 = 0;
                for (var i = 0; i < _envelope.Length; i++)
                {
                    if (_rankingMatrix[cell.Row][cell.Col][i] == 0) continue;

                    r1 += _rankingMatrix[cell.Row][cell.Col][i];
                    n1++;
                }

                var p1 = FitScoreCalculator.GetRankSumPvalue(_cachedNumPeaks[cell.Col], n1, r1);
                if (p1 < minPvalue) minPvalue = p1;
            }

            if (minPvalue < 10E-14) return 30;
            return -Math.Log(minPvalue);
        }

        private double CalculateHyperGeometricScore(ChargeLcScanCluster cluster)
        {
            if (_cachedXic == null) return 0;

            var n = (_cachedMaxBinNum - _cachedMinBinNum + 1) * cluster.Members.Count;
            var k = cluster.Members.Sum(cell => _cachedNumPeaks[cell.Col]);
            var n1 = cluster.Members.Count*_envelope.Length;
            var k1 = 0;
            foreach (var cell in cluster.Members)
            {
                for (var i = 0; i < _envelope.Length; i++)
                {
                    if (_xicMatrix[cell.Row][cell.Col][i] > 0) k1++;
                }
            }
            var pValue = FitScoreCalculator.GetHyperGeometricPvalue(n, k, n1, k1);

            if (pValue < 10E-200) return 450;
            return -Math.Log(pValue);
        }

        private double CalculateXicCorrelationOverTimeBetweenIsotopes(ChargeLcScanCluster cluster)
        {
            double ret = 0;
            var maxCol = cluster.MaxCol;
            var minCol = cluster.MinCol;
            var maxRow = cluster.MaxRow;
            var minRow = cluster.MinRow;
            
            var colLen = maxCol - minCol + 1;
            const int minColLen = 30;

            if (colLen < minColLen)
            {
                minCol = Math.Max(minCol - (int)((minColLen - colLen) * 0.5), 0);

                if (minCol == 0) maxCol = minCol + minColLen - 1;
                else
                {
                    maxCol = Math.Min(maxCol + (int)((minColLen - colLen) * 0.5), _nScans - 1);
                    if (maxCol == _nScans - 1)  minCol = maxCol - minColLen + 1;
                }
                colLen = maxCol - minCol + 1;
            }

            var xicProfile = new double[_topEnvelopes.Length][];
            for (var i = 0; i < xicProfile.Length; i++) xicProfile[i] = new double[colLen];

            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    if (_intensityMapFull[row][col] < 1E-6) continue;

                    for (var k = 0; k < _topEnvelopes.Length; k++)
                    {
                        xicProfile[k][col - minCol] += _xicMatrix[row][col][_topEnvelopes[k]];    
                    }
                }
            }

            for (var i = 1; i < xicProfile.Length; i++)
            {
                for (var j = i + 1; j < xicProfile.Length; j++)
                {
                    var coeff = FitScoreCalculator.GetPearsonCorrelation(xicProfile[i], xicProfile[j]);
                    if (coeff < 1.0 && coeff > ret)
                    {
                        ret = coeff;
                    }
                }
            }
            return ret;
        }

        private double CalculateXicCorrelationOverTimeBetweenCharges(ChargeLcScanCluster cluster)
        {
            var maxCol = cluster.MaxCol;
            var minCol = cluster.MinCol;
            var maxRow = cluster.MaxRow;
            var minRow = cluster.MinRow;

            var colLen = maxCol - minCol + 1;
            const int minColLen = 20;

            if (colLen < minColLen)
            {
                minCol = Math.Max(minCol - (int)((minColLen - colLen) * 0.5), 0);

                if (minCol == 0) maxCol = minCol + minColLen - 1;
                else
                {
                    maxCol = Math.Min(maxCol + (int)((minColLen - colLen) * 0.5), _nScans - 1);
                    if (maxCol == _nScans - 1) minCol = maxCol - minColLen + 1;
                }
                colLen = maxCol - minCol + 1;
            }
          
            if (Math.Abs(minRow - maxRow) < 5)
            {
                minRow = Math.Max(minRow - 3, 0);
                maxRow = Math.Min(maxRow + 3, _nCharges-1);
            }

            double ret = 0;
            for (var row1 = minRow; row1 <= maxRow; row1++)
            {
                for (var row2 = row1+1; row2 <= maxRow; row2++)
                {
                    var c = FitScoreCalculator.GetPearsonCorrelation(_intensityMapMostAbundant[row1], minCol, _intensityMapMostAbundant[row2], minCol, maxCol - minCol + 1);

                    if (c > ret) ret = c;
                }
            }
            return ret;
        }

        private int[] _topEnvelopes;

        private void SetIsoEnvelope()
        {
            const double envelopeTh = 0.09;
            var isoEnv = Averagine.GetIsotopomerEnvelope(_monoIsotopicMass);
            var isotopeRankings = ArrayUtil.GetRankings(isoEnv.Envolope);

            var nEnvelope = Math.Min(isoEnv.Envolope.Count(x => x > envelopeTh), MaxEnvelopeLength);
            var nTopEnvelope = Math.Min(nEnvelope, 4);

            _envelope           = new double[nEnvelope];
            _envelopeIndex      = new int[nEnvelope];
            _topEnvelopes       = new int[nTopEnvelope];
            
            var index = 0;
            var indexTop = 0;
            for (var i = 0; i < isoEnv.Envolope.Length; i++)
            {
                if (isotopeRankings[i] > nEnvelope) continue;

                _envelope[index] = isoEnv.Envolope[i];
                _envelopeIndex[index] = i;
                
                if (i == isoEnv.MostAbundantIsotopeIndex) _mostAbundantIsotopeIndex = index;

                if (isotopeRankings[i] <= nTopEnvelope)
                {
                    _topEnvelopes[indexTop] = index;
                    indexTop++;
                }

                index++;
            }
        }

        private void BuildMatrix()
        {
            //var tolerance = new Tolerance(20);
            SetIsoEnvelope();

            Array.Clear(_observedCharges, 0, _observedCharges.Length);

            for (var c = 0; c < _nCharges; c++)
            {
                Array.Clear(_intensityMapFull[c], 0, _nScans);
                Array.Clear(_intensityMapMostAbundant[c], 0, _nScans);
                Array.Clear(_correlationMap[c], 0, _nScans);

                //for (var j = 0; j < _nScans; j++) Array.Clear(_xicMatrix[c][j], 0, MaxEnvelopeLength);
            }

            // Fill _intensityMap and _xicMatrix arrays
            for (var c = 0; c < _nCharges; c++)
            {
                for (var i = 0; i < _envelopeIndex.Length; i++)
                {
                    var isotopeIndex = _envelopeIndex[i];
                    var curMz = Ion.GetIsotopeMz(_monoIsotopicMass, c + _minCharge, isotopeIndex);
                    
                    var curBinNum = _comparer.GetBinNumber(curMz);
                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum)
                    {
                        for (var j = 0; j < _nScans; j++)
                        {
                            _xicMatrix[c][j][i] = 0;
                            _rankingMatrix[c][j][i] = 0;
                        }
                        continue;
                    }
                    
                    double[] xic = null;

                    if (_cachedXic != null)
                    {
                        xic = _cachedXic[curBinNum - _cachedMinBinNum];
                    }
                    else
                    {
                        var mzStart = _comparer.GetMzStart(curBinNum);
                        var mzEnd = _comparer.GetMzEnd(curBinNum);
                        var mzHalf = 0.5 * (mzEnd - mzStart);
                        xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart - mzHalf, mzEnd + mzHalf);

                        //xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(curMz - tolerance.GetToleranceAsTh(curMz), curMz + tolerance.GetToleranceAsTh(curMz));
                    }

                    if (!_observedCharges[c] && xic.Any(x => x > 0)) _observedCharges[c] = true;

                    for (var j = 0; j < _nScans; j++)
                    {
                        if (i == _mostAbundantIsotopeIndex) _intensityMapMostAbundant[c][j] += xic[j];
                        
                        _intensityMapFull[c][j] += xic[j];
                        _xicMatrix[c][j][i] = xic[j];

                        if (_cachedRanking != null)
                        {
                            _rankingMatrix[c][j][i] = _cachedRanking[j][curBinNum - _cachedMinBinNum];
                        }
                    }
                    /*
                    if (c < 30) continue;

                    var cmax = Math.Min(c + 3, _nCharges - 1);
                    var cmin = cmax - 6;
                    
                    for (var c2 = cmin; c2 <= cmax; c2++)
                    {
                        if (c == c2) continue;

                        for (var j = 0; j < _nScans; j++)
                        {
                            _intensityMap[c2][j] += xic[j];
                            _xicMatrix[c2][j][i] += xic[j];
                        }
                    }*/

                }
            }

            var k = 0;
            _chargeIndexes = new int[_observedCharges.Where((x) => x == true).Count()];
            for (var i = 0; i < _nCharges; i++)
            {
                if (!_observedCharges[i]) continue;
                _chargeIndexes[k] = i;
                k++;
            }

            // Fill _correlationMap array
            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (_intensityMapFull[i][j] > 0)
                    {
                        _correlationMap[i][j] = FitScoreCalculator.GetPearsonCorrelation(_xicMatrix[i][j], 0, _envelope, 0, _envelope.Length);
                    }
                    else
                    {
                        _correlationMap[i][j] = 0;
                    }
                }
            }
        }

        public double[][] GetIntensityMap()
        {
            return _intensityMapFull;
        }

        public double[][] GetCorrelationMap()
        {
            return _correlationMap;
        }
        
        private static double GetThreshold(double topRatio, IEnumerable<double[]> data)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            foreach (var t in data)
                nonZeroIntensities.AddRange(t.Where(t1 => t1 > 0));

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[(int)Math.Floor(nonZeroIntensities.Count * (1-topRatio))];

            return intensityThreshold;
        }

        private const double CorrSeedLowerBound = 0.5d;
        private const double CorrLowerBound = 0.5d;
        private IEnumerable<ChargeLcScanCell> GetSeedCells()
        {
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();

            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (_correlationMap[i][j] < CorrSeedLowerBound) continue;

                    seedCells.Add(new KeyValuePair<double, ChargeLcScanCell>(_correlationMap[i][j], new ChargeLcScanCell(i, j)));
                }
            }

            return seedCells.OrderByDescending(x => x.Key).Select(x => x.Value);
        }

        // cutoff < Corr(thoeretical_envelope) * 2 + Corr(between isotopes)
        private List<ChargeLcScanCluster> FindClusters(double cutoff = 2.5, double hypergeometricCutoff = 50, double ranksumCutoff = 4)
        {
            BuildMatrix();

            var checkedOut = new bool[_nCharges][];
            for (var i = 0; i < _nCharges; i++) checkedOut[i] = new bool[_nScans];
            
            var clusters = new List<ChargeLcScanCluster>();
            var tempEnvelope = new double[_envelope.Length];

            foreach (var seedCell in GetSeedCells())
            {
                if (checkedOut[seedCell.Row][seedCell.Col]) continue;

                var newCluster = new ChargeLcScanCluster(seedCell, _xicMatrix[seedCell.Row][seedCell.Col], _envelope.Length, _intensityMapFull[seedCell.Row][seedCell.Col], _correlationMap[seedCell.Row][seedCell.Col]);
                var neighbors   = new Queue<ChargeLcScanCell>();

                neighbors.Enqueue(seedCell); // pick a seed
                checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    const int chargeNeighborGap = 2;
                    
                    for (var k = Math.Max(cell.Row - chargeNeighborGap, _chargeIndexes.First()); k <= Math.Min(cell.Row + chargeNeighborGap, _chargeIndexes.Last()); k++)
                    {
                        for (var l = Math.Max(cell.Col - 1, 0); l <= Math.Min(cell.Col + 1, _nScans-1); l++)
                        {
                            if (checkedOut[k][l] || _correlationMap[k][l] < CorrLowerBound) continue;
                            
                            for (var t = 0; t < tempEnvelope.Length; t++) tempEnvelope[t] = _xicMatrix[cell.Row][cell.Col][t] + _xicMatrix[k][l][t];

                            var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelope, tempEnvelope);

                            if (_correlationMap[cell.Row][cell.Col] < newScore || _correlationMap[k][l] > 0.85)
                            {
                                var newMember = new ChargeLcScanCell(k, l);
                                neighbors.Enqueue(newMember);
                                newCluster.AddMember(newMember, _intensityMapFull[k][l], _xicMatrix[k][l], newScore);
                                checkedOut[k][l] = true;
                            }
                        }
                    }
                }

                newCluster.CorrBetweenObsTh = FitScoreCalculator.GetPearsonCorrelation(_envelope, newCluster.ObservedEnvelope);

                newCluster.HyperGeometricScore = CalculateHyperGeometricScore(newCluster);
                if (newCluster.HyperGeometricScore < hypergeometricCutoff) continue;

                newCluster.RankSumScore = CalculateRankSumScore(newCluster);
                if (newCluster.RankSumScore - 0.0005 * _monoIsotopicMass < ranksumCutoff) continue;
                
                newCluster.CorrBetweenIsotopes = CalculateXicCorrelationOverTimeBetweenIsotopes(newCluster);
                if (newCluster.GetScore() < cutoff) continue;

                clusters.Add(newCluster);

                for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) checkedOut[i][j] = true;
                
            }

            clusters.Sort(ChargeLcScanCluster.CompareByCorrelation);
            return clusters;
        }

        /*
        public void PreparePrecursorXicVector()
        {
            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);
            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];

            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                var mzStart = _comparer.GetMzStart(binNum);
                var mzEnd = _comparer.GetMzEnd(binNum);
                var mzHalf = 0.5*(mzEnd - mzStart);
                var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart - mzHalf, mzEnd + mzHalf);

                _cachedXic[binNum - _cachedMinBinNum] = xic;
            }
        }
                
        public void PrepareSmoothFullPrecursorXicVector()
        {
            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);
            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];

            var binNumToAvgMass = new double[_cachedMaxBinNum - _cachedMinBinNum + 1];
            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                _cachedXic[binNum - _cachedMinBinNum] = new double[_nScans];
                binNumToAvgMass[binNum - _cachedMinBinNum] = _comparer.GetMzAverage(binNum);
            }

            var xic = ((PbfLcMsRun)_run).GetAllPrecursorExtractedIonChromatogram();
            var tolerance = _comparer.GetTolerance();

            foreach (var xicPoint in xic)
            {
                var xicBinNum = _comparer.GetBinNumber(xicPoint.Mz);

                var bandWidth = tolerance.GetToleranceAsTh(xicPoint.Mz) * 0.5;
                var minBinNum = Math.Max(xicBinNum - 1, _cachedMinBinNum);
                var maxBinNum = Math.Min(xicBinNum + 1, _cachedMaxBinNum);
               
                for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                {
                    var u = (xicPoint.Mz - binNumToAvgMass[binNum - _cachedMinBinNum]) / bandWidth;
                    var smoothedIntensity = xicPoint.Intensity * (Math.Exp(-Math.Pow(u, 2) / 2));
                    _cachedXic[binNum - _cachedMinBinNum][_ms1ScanNumToIndex[xicPoint.ScanNum]] += smoothedIntensity;
                }
            }
        }
        */
   }

    internal class ChargeLcScanCluster
    {
        //internal double Score { get; private set; }

        internal double CorrBetweenObsTh;
        internal double CorrBetweenIsotopes;
        internal double RankSumScore;
        internal double HyperGeometricScore;

        internal readonly List<ChargeLcScanCell> Members;

        internal double HighestIntensity { get; private set; }
        internal double TotalIntensity { get; private set; }

        private int _highestIntensityMemberIndex;
        internal double[] ObservedEnvelope { get; private set; }
        internal double SeedCorrelation { get; private set; }

        internal int MaxCol { get; private set; }
        internal int MinCol { get; private set; }
        internal int MaxRow { get; private set; }
        internal int MinRow { get; private set; }

        internal bool Active;
        internal double[] IntensityAlongCol;

        internal ChargeLcScanCluster(ChargeLcScanCell seed, double[] seedEnvelope, int envelopeLength, double seedIntensity, double seedCorr)
        {
            CorrBetweenObsTh = seedCorr;
            CorrBetweenIsotopes = 0;
            RankSumScore = 0;
            HyperGeometricScore = 0;

            SeedCorrelation = seedCorr;
            Members = new List<ChargeLcScanCell> {seed};

            ObservedEnvelope = new double[envelopeLength];
            Array.Copy(seedEnvelope, ObservedEnvelope, envelopeLength);

            HighestIntensity = seedIntensity;
            _highestIntensityMemberIndex = 0;

            MaxCol = seed.Col;
            MinCol = seed.Col;
            MaxRow = seed.Row;
            MinRow = seed.Row;
            TotalIntensity = seedIntensity;
        }

        internal ChargeScanRange GetChargeScanRange(int[] ms1ScanNums, int minCharge)
        {
            var seedScanNum = ms1ScanNums[Members[0].Col];
            return new ChargeScanRange(MinRow + minCharge, MaxRow + minCharge, ms1ScanNums[MinCol], ms1ScanNums[MaxCol], seedScanNum, GetScore());
        }
        /*
        internal void SetScore(double s)
        {
            Score = s;
        }*/

        public double GetScore()
        {
            return CorrBetweenObsTh*3 + CorrBetweenIsotopes;
        }

        public bool Overlaps(ChargeLcScanCluster other)
        {
            return ((MinCol <= other.MinCol && other.MinCol <= MaxCol) ||
                    (MinCol <= other.MaxCol && other.MaxCol <= MaxCol)) &&
                   ((MinRow <= other.MinRow && other.MinRow <= MaxRow) ||
                    (MinRow <= other.MaxRow && other.MaxRow <= MaxRow));
        }

        internal ChargeLcScanCell GetHighestIntensityCell()
        {
            return Members[_highestIntensityMemberIndex];
        }

        internal void AddMember(ChargeLcScanCell member, double memberIntensity, double[] newObservedEnvelope, double newCorr)
        {
            if (memberIntensity > HighestIntensity)
            {
                HighestIntensity = memberIntensity;
                _highestIntensityMemberIndex = Members.Count;
            }

            TotalIntensity += memberIntensity;

            if (member.Col > MaxCol) MaxCol = member.Col;
            else if (member.Col < MinCol) MinCol = member.Col;
            if (member.Row > MaxRow) MaxRow = member.Row;
            else if (member.Row < MinRow) MinRow = member.Row;

            CorrBetweenObsTh = newCorr;
            for (var i = 0; i < ObservedEnvelope.Length; i++) ObservedEnvelope[i] += newObservedEnvelope[i];

            Members.Add(member);
        }

        internal static int CompareByCorrelation(ChargeLcScanCluster x, ChargeLcScanCluster y)
        {
            return -x.GetScore().CompareTo(y.GetScore());
        }
    }

    class ChargeLcScanCell
    {
        internal int Row;
        internal int Col;

        internal ChargeLcScanCell(int row, int col)
        {
            Row = row;
            Col = col;
        }
    }
    
}
