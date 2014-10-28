using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.Remoting.Messaging;
using System.Text;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using MultiDimensionalPeakFinding;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {
        public ChargeLcScanMatrix(LcMsRun run, int numBits = 27, bool useCache = true)
        {
            _minCharge = 2;
            _maxCharge = 50;
            _ms1ScanNums = run.GetMs1ScanVector();
            _comparer = new MzComparerWithBinning(numBits);
            _smoother = new SavitzkyGolaySmoother(9, 2);

            _nScans = _ms1ScanNums.Length;
            _nCharges = _maxCharge - _minCharge + 1;

            _run = run;
            _clusterMap = new Dictionary<int, IEnumerable<ChargeScanRange>>();
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
            }

            _intensityMap       = new double[_nCharges][];
            _smoothIntensityMap = new double[_nCharges][];
            _correlationMap     = new double[_nCharges][];
            _xicMatrix          = new double[_nCharges][][];
            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMap[i] = new double[_nScans];
                _smoothIntensityMap[i] = new double[_nScans];
                _correlationMap[i] = new double[_nScans];
                _xicMatrix[i] = new double[_nScans][];
                for (var j = 0; j < _nScans; j++)
                {
                    _xicMatrix[i][j] = new double[MaxEnvelopeLength];
                }
            }
        }

        public IEnumerable<ChargeScanRange> GetProbableChargeScanRegions(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);

            IEnumerable<ChargeScanRange> ranges;
            if (_clusterMap.TryGetValue(binNumber, out ranges)) return ranges;
            
            BuildMatrix(binNumber);
            var clusters = FindClusters(0.7, 0.6);
            
            ranges = clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
            return _clusterMap[binNumber] = ranges;
            
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double monoIsotopicMass)
        {
            IEnumerable<int> ms2Scans;
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            if (_binNumToMs2ScanNumsMap.TryGetValue(binNumber, out ms2Scans)) return ms2Scans;

            var ms2ScanList = new List<int>();
            var isoEnv = Averagine.GetIsotopomerEnvelope(monoIsotopicMass);
            var mostAbundantIsotopeIndex = isoEnv.MostAbundantIsotopeIndex;

            foreach (var chargeScanRange in GetProbableChargeScanRegions(monoIsotopicMass))
            {
                var minScanNum = chargeScanRange.MinScanNum;
                var maxScanNum = chargeScanRange.MaxScanNum;
                for (var charge = chargeScanRange.MinCharge; charge <= chargeScanRange.MaxCharge; charge++)
                {
                    var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoIsotopicMass, charge, mostAbundantIsotopeIndex);
                    ms2ScanList.AddRange(_run.GetFragmentationSpectraScanNums(mostAbundantIsotopeMz)
                        .Where(scanNum => scanNum >= minScanNum && scanNum <= maxScanNum));
                }
            }
            ms2ScanList.Sort();

            _binNumToMs2ScanNumsMap.Add(binNumber, ms2ScanList);
            return ms2ScanList;
        }

        public void SetChargeRange(int minCharge, int maxCharge)
        {
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _nCharges = _maxCharge - _minCharge + 1;
        }

        internal IEnumerable<ChargeLcScanCluster> GetProbableChargeScanCluster(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
           
            BuildMatrix(binNumber);
            var clusters = FindClusters(0.7, 0.6);

            foreach (var c in clusters)
            {
                c.IntensityAlongCol = new double[c.MaxCol - c.MinCol + 1];
                for (var col = c.MinCol; col <= c.MaxCol; col++)
                {
                    for (var row = c.MinRow; row <= c.MaxRow; row++)
                    {
                        c.IntensityAlongCol[col-c.MinCol] += _intensityMap[row][col];
                    }
                }
            }

            return clusters;
        }

        public IEnumerable<ChargeScanRange> GetAllChargeScanRegions(double proteinMass, out double[][] scores)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);
            //BuildMatrix(proteinMass);
            var clusters = FindClusters(0.01, 0.01);

            scores = new double[4][];
            for (var i = 0; i < scores.Length; i++) scores[i] = new double[clusters.Count];
            
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            var mostAbundantIsotopeIndex = isoEnv.MostAbundantIsotopeIndex;
            
            for(var i = 0; i < clusters.Count; i++)
            {
                scores[0][i] = clusters[i].Score; //FitScoreCalculator.GetPearsonCorrelation(_envelope, 0, clusters[i].ObservedEnvelope, 0, _envelope.Length);
                scores[1][i] = CalculateXicCorrelationOverTimeBetweenCharges(clusters[i]);
                scores[2][i] = CalculateXicCorrelationOverTimeBetweenIsotopes(proteinMass, clusters[i]);
                var chargeScanRange = clusters[i].GetChargeScanRange(_ms1ScanNums, _minCharge);
                var minScanNum = chargeScanRange.MinScanNum;
                var maxScanNum = chargeScanRange.MaxScanNum;
                for (var charge = chargeScanRange.MinCharge; charge <= chargeScanRange.MaxCharge; charge++)
                {
                    var mostAbundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, mostAbundantIsotopeIndex);
                    scores[3][i] +=  _run.GetFragmentationSpectraScanNums(mostAbundantIsotopeMz).Where(scanNum => scanNum >= minScanNum && scanNum <= maxScanNum).Count();
                }
            }
            
            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }

        public MzComparerWithBinning GetMzComparerWithBinning()
        {
            return _comparer;
        }

        private readonly LcMsRun _run;
        private int _minCharge;
        private int _maxCharge;
        private int _nCharges;
        private readonly int _nScans;
        private readonly int[] _ms1ScanNums;
        private readonly double[][] _intensityMap;
        private readonly double[][] _correlationMap;
        private readonly double[][][] _xicMatrix;

        private readonly double[][] _smoothIntensityMap;
        private double[] _envelope;
        private int[] _envelopeIndex = null;
        private const int MaxEnvelopeLength = 50;
        private int[] _chargeIndexes;
        
        private readonly MzComparerWithBinning _comparer;
        private readonly SavitzkyGolaySmoother _smoother;

        // caching
        private readonly double[][] _cachedXic;
        private readonly int _cachedMinBinNum;
        private readonly int _cachedMaxBinNum;

        private readonly Dictionary<int, IEnumerable<ChargeScanRange>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;
        
        private double CalculateXicCorrelationOverTimeBetweenIsotopes(double proteinMass, ChargeLcScanCluster cluster)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            double ret = 0;
            var maxCol = cluster.MaxCol;
            var minCol = cluster.MinCol;
            //var maxRow = cluster.MaxRow;
            //var minRow = cluster.MinRow;

            var minRow = _chargeIndexes.First();
            var maxRow = _chargeIndexes.Last();
            
            var colLen = maxCol - minCol + 1;
            //int minColLen = Math.Max((int)Math.Ceiling(colLen * 1.5), 30);
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

            var mostAbundantIsotopeIndex = 0;
            /*
            double maxEnv = 0;
            for (var i = 0; i < _envelope.Length; i++)
            {
                if (_envelope[i] > maxEnv)
                {
                    maxEnv = _envelope[i];
                    mostAbundantIsotopeIndex = i;
                }
            }*/
            const int isotopeDeltaIndex = 1;
            
            for (var i = 0; i < _envelopeIndex.Length; i++)
            {
                if (_envelopeIndex[i] != isoEnv.MostAbundantIsotopeIndex) continue;
                mostAbundantIsotopeIndex = i;
                break;
            }
               /*         
             var charge = minRow + _minCharge;
             var mostAbundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, mostAbundantIsotopeIndex);

             var binNum = _comparer.GetBinNumber(mostAbundantIsotopeMz);
             var binMzSize = _comparer.GetMzEnd(binNum) - _comparer.GetMzStart(binNum);

             var minIsotopeDeltaIndex = (charge*binMzSize*1.1)/Biology.Constants.C13MinusC12;
             var isotopeDeltaIndex = (int) Math.Ceiling(minIsotopeDeltaIndex);
             */
            var xicProfile = new double[4][];
            for (var i = 0; i < xicProfile.Length; i++) xicProfile[i] = new double[colLen];

            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    if (_intensityMap[row][col] < 1E-6) continue;

                    xicProfile[0][col - minCol] += _xicMatrix[row][col][mostAbundantIsotopeIndex];
                    xicProfile[1][col - minCol] += _xicMatrix[row][col][mostAbundantIsotopeIndex - isotopeDeltaIndex];
                    xicProfile[2][col - minCol] += _xicMatrix[row][col][mostAbundantIsotopeIndex + isotopeDeltaIndex];
                    xicProfile[3][col - minCol] += _xicMatrix[row][col][mostAbundantIsotopeIndex + isotopeDeltaIndex*2];
                }
            }

            for (var i = 0; i < xicProfile.Length; i++)
                xicProfile[i] = _smoother.Smooth(xicProfile[i]);

            /*
            var temp = new List<KeyValuePair<double, int>>();
            for (var i = 0; i < xicProfile.Length; i++)
            {
                xicProfile[i] = _smoother.Smooth(xicProfile[i]);

                temp.Add(new KeyValuePair<double, int>(xicProfile[i].Sum(), i));
            }

            temp = temp.OrderByDescending(x => x.Key).ToList();
            ret = Math.Max(ret, FitScoreCalculator.GetPearsonCorrelation(xicProfile[temp[0].Value], xicProfile[temp[1].Value]));*/


            for (var i = 0; i < xicProfile.Length; i++)
            {
                for (var j = i + 1; j < xicProfile.Length; j++)
                {
                    ret = Math.Max(ret, FitScoreCalculator.GetPearsonCorrelation(xicProfile[i], xicProfile[j]));
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
            //int minColLen = Math.Max((int)Math.Ceiling(colLen * 1.5), 30);
            const int minColLen = 30;

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
            /*
            var row1 = 0;
            var row2 = 0;
            if (maxRow == minRow)
            {
                row1 = minRow;
                row2 = (row1 == _intensityMap.Length - 1) ? row1 - 1 : row1 + 1;
            }
            else
            {
                var temp = new List<KeyValuePair<double, int>>();
                for (var row = minRow; row <= maxRow; row++)
                {
                    var intensity = 0d;
                    for (var col = minCol; col <= maxCol; col++)
                    {
                        intensity += _intensityMap[row][col];
                    }
                    temp.Add(new KeyValuePair<double, int>(intensity, row));
                }

                var i = 0;
                foreach (var t in temp.OrderByDescending(x => x.Key))
                {
                    if (i == 0) row1 = t.Value;
                    else
                    {
                        row2 = t.Value;
                        break;
                    }
                    i++;
                }
            }

            return FitScoreCalculator.GetPearsonCorrelation(_smoothIntensityMap[row1], minCol, _smoothIntensityMap[row2], minCol, colLen);
            */

            var row1 = (int)Math.Floor((maxRow + minRow) * 0.5);
            var row2 = row1 + 1;
            var row3 = row1 - 1;
            if (row2 == _intensityMap.Length)
            {
                row1--;
                row2--;
                row3--;
            }
            else if (row3 < 0)
            {
                row1++;
                row2++;
                row3++;
            }

            return 0.5*
                   (FitScoreCalculator.GetPearsonCorrelation(_intensityMap[row1], minCol, _intensityMap[row2], minCol,
                       maxCol - minCol + 1) +
                    FitScoreCalculator.GetPearsonCorrelation(_intensityMap[row1], minCol, _intensityMap[row3], minCol,
                        maxCol - minCol + 1));
        }

        private void BuildMatrix(int binNum)
        {
            double proteinMass = _comparer.GetMzAverage(binNum);
            BuildMatrix(proteinMass);
        }
        
        private void BuildMatrix(double proteinMass)
        {
            const double epsillon = 1E-6;

            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            var n = isoEnv.Envolope.Count(x => x > 0.1);
            if (n <= MaxEnvelopeLength)
            {
                _envelopeIndex = new int[n];
                _envelope = new double[n];
                var idx = 0;
                for (var i = 0; i < isoEnv.Envolope.Length; i++)
                {
                    if (isoEnv.Envolope[i] > 0.1)
                    {
                        _envelopeIndex[idx] = i;
                        _envelope[idx] = isoEnv.Envolope[i];
                        idx++;
                    }
                }
            }
            else
            {
                _envelope = new double[MaxEnvelopeLength];
                _envelopeIndex = new int[MaxEnvelopeLength];
                var sortedEnvelope = isoEnv.Envolope.ToList().OrderByDescending(x => x).ToArray();
                var i = 0;
                for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                {
                    if (isoEnv.Envolope[isotopeIndex] > sortedEnvelope[MaxEnvelopeLength])
                    {
                        _envelopeIndex[i] = isotopeIndex;
                        _envelope[i] = isoEnv.Envolope[isotopeIndex];
                        i++;
                    }
                }
            }

            var observedCharges = new bool[_nCharges];

            // Fill _intensityMap and _xicMatrix arrays
            for (var charge = _minCharge; charge <= _maxCharge; charge++)
            {
                Array.Clear(_intensityMap[charge - _minCharge], 0, _nScans);
                //Array.Clear(_correlationMap[charge - _minCharge], 0, _nScans);

                for (var i = 0; i < _envelopeIndex.Length; i++)
                {
                    var isotopeIndex = _envelopeIndex[i];
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);
                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum)
                    {
                        for (var j = 0; j < _nScans; j++)
                        {
                            _xicMatrix[charge - _minCharge][j][i] = 0;
                        }
                        continue;
                    }
                    double[] xic = null;

                    if (_cachedXic != null)
                        xic = _cachedXic[curBinNum - _cachedMinBinNum];
                    else
                    {
                        var mzStart = _comparer.GetMzStart(curBinNum);
                        var mzEnd = _comparer.GetMzEnd(curBinNum);
                        var mzHalf = 0.5 * (mzEnd - mzStart);
                        xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart - mzHalf, mzEnd + mzHalf);
                    }

                    if (!observedCharges[charge - _minCharge] && xic.Any(x => x > 0)) observedCharges[charge - _minCharge] = true;

                    for (var j = 0; j < _nScans; j++)
                    {
                        _intensityMap[charge - _minCharge][j] += xic[j];
                        _xicMatrix[charge - _minCharge][j][i] = xic[j];
                    }
                }
            }

            var k = 0;
            _chargeIndexes = new int[observedCharges.Where((x) => x == true).Count()];
            for (var i = 0; i < _nCharges; i++)
            {
                if (observedCharges[i])
                {
                    _chargeIndexes[k] = i;
                    k++;
                }
            }

            // Fill _correlationMap array
            foreach (var i in _chargeIndexes)
            {
                _smoothIntensityMap[i] = _smoother.Smooth(_intensityMap[i]);
                
                for (var j = 0; j < _nScans; j++)
                {
                    if (_intensityMap[i][j] < epsillon)
                    {
                        _correlationMap[i][j] = 0;
                    }
                    else
                    {
                        _correlationMap[i][j] = FitScoreCalculator.GetPearsonCorrelation(_xicMatrix[i][j], 0, _envelope, 0,
                                                                                        _envelope.Length);
                    }
                }
            }
        }

        public double[][] GetIntensityMap()
        {
            return _intensityMap;
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
        
        private IEnumerable<ChargeLcScanCell> GetSeedCells(double corrLowerBound = 0.5)
        {
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();

            foreach (var i in _chargeIndexes)
            {
                for (var j = 0; j < _nScans; j++)
                {
                    if (_correlationMap[i][j] < corrLowerBound) continue;

                    seedCells.Add(new KeyValuePair<double, ChargeLcScanCell>(_correlationMap[i][j], new ChargeLcScanCell(i, j)));
                }
            }

            return seedCells.OrderByDescending(x => x.Key).Select(x => x.Value);
        }

        private List<ChargeLcScanCluster> FindClusters(double envelopCorrTh = 0.7, double chargeCorrTh = 0.5)
        {
            //const double corrLowerBound = 0.1d;
            //var corrSeedLowerBound = GetThreshold(0.5, _correlationMap); // median value in correlations
            //var corrLowerBound = corrSeedLowerBound;
            const double corrSeedLowerBound = 0.5d;
            const double corrLowerBound = 0.5d;


            var checkedOut = new bool[_nCharges][];
            for (var i = 0; i < _nCharges; i++) checkedOut[i] = new bool[_nScans];
            
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seedCell in GetSeedCells(corrSeedLowerBound))
            {
                if (checkedOut[seedCell.Row][seedCell.Col]) continue;
                var newCluster = new ChargeLcScanCluster(seedCell, _xicMatrix[seedCell.Row][seedCell.Col], _intensityMap[seedCell.Row][seedCell.Col], _correlationMap[seedCell.Row][seedCell.Col]);

                var neighbors = new Queue<ChargeLcScanCell>();
                neighbors.Enqueue(seedCell); // pick a seed
                checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    //var chargeNeighborRange = (int) Math.Ceiling(cell.Row*0.1);
                    const int chargeNeighborRange = 2;

                    for (var k = Math.Max(cell.Row - chargeNeighborRange, _chargeIndexes.First()); k <= Math.Min(cell.Row + chargeNeighborRange, _chargeIndexes.Last()); k++)
                    {
                        for (var l = Math.Max(cell.Col - 2, 0); l <= Math.Min(cell.Col + 2, _nScans-1); l++)
                        {
                            if (checkedOut[k][l] || _intensityMap[k][l] < 10E-6 || _correlationMap[k][l] < corrLowerBound) continue;

                            var tempEnvelope = new double[_envelope.Length];
                            for (var t = 0; t < tempEnvelope.Length; t++) tempEnvelope[t] = newCluster.ObservedEnvelope[t] + _xicMatrix[k][l][t];

                            var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelope, tempEnvelope);

                            if (newScore >= newCluster.Score || _correlationMap[k][l] > 0.9)
                            {
                                var newMember = new ChargeLcScanCell(k, l);
                                neighbors.Enqueue(newMember);
                                newCluster.AddMember(newMember, _intensityMap[k][l], tempEnvelope, newScore);
                                checkedOut[k][l] = true;

                                /*for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                                    for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) checkedOut[i][j] = true;*/
                            }
                        }
                    }
                }

                if (newCluster.Score > envelopCorrTh)
                {
                    var score2 = CalculateXicCorrelationOverTimeBetweenCharges(newCluster);
                    if (score2 > chargeCorrTh || newCluster.Score > 0.88) clusters.Add(newCluster);

                    for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                        for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) checkedOut[i][j] = true;
                }
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
        
        public void PrepareFullPrecursorXicVector2()
        {
            
            var _ms1ScanNumToIndex = new int[_ms1ScanNums[_ms1ScanNums.Length - 1] + 1];
            for (var i = 0; i < _ms1ScanNums.Length; i++) _ms1ScanNumToIndex[_ms1ScanNums[i]] = i;

            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);
            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];
            
            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                _cachedXic[binNum - _cachedMinBinNum] = new double[_nScans];
            }

            var xic = ((PbfLcMsRun)_run).GetAllPrecursorExtractedIonChromatogram();

            foreach (var xicPoint in xic)
            {
                var xicBinNum = _comparer.GetBinNumber(xicPoint.Mz);

                var minBinNum = 0;
                var maxBinNum = 0;
                if (xicPoint.Mz < _comparer.GetMzAverage(xicBinNum))
                {
                    minBinNum = Math.Max(xicBinNum - 1, _cachedMinBinNum);
                    maxBinNum = xicBinNum;
                }
                else
                {
                    minBinNum = xicBinNum;
                    maxBinNum = Math.Min(xicBinNum + 1, _cachedMaxBinNum);                    
                }

                for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                {
                    _cachedXic[binNum - _cachedMinBinNum][_ms1ScanNumToIndex[xicPoint.ScanNum]] = Math.Max(_cachedXic[binNum - _cachedMinBinNum][_ms1ScanNumToIndex[xicPoint.ScanNum]], xicPoint.Intensity);
                }
            }
        }
        */
   }

    internal class ChargeLcScanCluster
    {
        internal double Score { get; private set; }
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

        internal ChargeLcScanCluster(ChargeLcScanCell seed, double[] seedEnvelope, double seedIntensity, double seedScore)
        {
            Score = seedScore;
            SeedCorrelation = seedScore;
            Members = new List<ChargeLcScanCell> {seed};
            ObservedEnvelope = seedEnvelope;
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
            return new ChargeScanRange(MinRow + minCharge, MaxRow + minCharge, ms1ScanNums[MinCol], ms1ScanNums[MaxCol], seedScanNum, Score);
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

        internal void AddMember(ChargeLcScanCell member, double memberIntensity, double[] newObservedEnvelope, double newScore)
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

            Score = newScore;
            ObservedEnvelope = newObservedEnvelope;
            Members.Add(member);
        }

        internal static int CompareByCorrelation(ChargeLcScanCluster x, ChargeLcScanCluster y)
        {
            return -x.Score.CompareTo(y.Score);
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
