using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {

        public ChargeLcScanMatrix(LcMsRun run, int numBits = 27)
        {
            _minCharge = 2;
            _maxCharge = 50;
            _ms1ScanNums = run.GetMs1ScanVector();
            _comparer = new MzComparerWithBinning(numBits);

            _nScans = _ms1ScanNums.Length;
            _nCharges = _maxCharge - _minCharge + 1;

            _run = run;
            _clusterMap = new Dictionary<int, IEnumerable<ChargeScanRange>>();
            _binNumToMs2ScanNumsMap = new Dictionary<int, IEnumerable<int>>();
            
            //Cache Xic data
            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);
            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];

            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                var mzStart = _comparer.GetMzStart(binNum);
                var mzEnd = _comparer.GetMzEnd(binNum);
                var mzHalf = 0.5 * (mzEnd - mzStart);
                var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart - mzHalf, mzEnd + mzHalf);
                _cachedXic[binNum - _cachedMinBinNum] = xic;
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
        

        public IEnumerable<ChargeScanRange> GetAllChargeScanRegions(double proteinMass, out double[][] scores)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);
            //BuildMatrix(proteinMass);

            var clusters = FindClusters(0.01, 0.01);
            clusters.Sort(ChargeLcScanCluster.CompareByCorrelation);

            scores = new double[4][];
            for (var i = 0; i < scores.Length; i++) scores[i] = new double[clusters.Count];
            
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            var mostAbundantIsotopeIndex = isoEnv.MostAbundantIsotopeIndex;
            
            for(var i = 0; i < clusters.Count; i++)
            {
                scores[0][i] = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, clusters[i].ObservedEnvelope);
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

        private readonly LcMsRun _run;
        private int _minCharge;
        private int _maxCharge;
        private int _nCharges;
        private readonly int _nScans;

        private double[][] _intensityMap;
        private double[][] _correlationMap;
        private double[][][] _xicMatrix;

        private int[] _chargeIndexes;

        private double[] _envelopePdf;
        private readonly int[] _ms1ScanNums;

        private readonly MzComparerWithBinning _comparer;

        // caching
        private readonly double[][] _cachedXic;
        private readonly int _cachedMinBinNum;
        private readonly int _cachedMaxBinNum;

        private readonly Dictionary<int, IEnumerable<ChargeScanRange>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;

        private double CalculateXicCorrelationOverTimeBetweenIsotopes(double proteinMass, ChargeLcScanCluster cluster)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            var maxCol = cluster.Members.Select(m => m.Col).Max();
            var minCol = cluster.Members.Select(m => m.Col).Min();
            var maxRow = cluster.Members.Select(m => m.Row).Max();
            var minRow = cluster.Members.Select(m => m.Row).Min();
            
            var colLen = maxCol - minCol + 1;
            const int minColLen = 20;

            if (colLen < minColLen)
            {
                minCol = Math.Max(minCol - (int)((minColLen - colLen) * 0.5), 0);

                if (minCol == 0) maxCol = minCol + minColLen - 1;
                else
                {
                    maxCol = Math.Min(maxCol + (int)((minColLen - colLen) * 0.5), _nScans - 1);
                    if (maxCol == _nScans - 1)  minCol = maxCol - minColLen + 1;
                }
            }

            var profile1 = new double[maxCol - minCol + 1];
            var profile2 = new double[maxCol - minCol + 1];
            var profile3 = new double[maxCol - minCol + 1];

            for (var row = minRow; row <= maxRow; row++)
            {
                for (var col = minCol; col <= maxCol; col++)
                {
                    if (_intensityMap[row][col] < 1) continue;

                    profile1[col - minCol] += _xicMatrix[row][col][isoEnv.MostAbundantIsotopeIndex];
                    profile2[col - minCol] += _xicMatrix[row][col][isoEnv.MostAbundantIsotopeIndex - 1];
                    profile3[col - minCol] += _xicMatrix[row][col][isoEnv.MostAbundantIsotopeIndex + 1];
                }
            }

            return 0.5*
                   (FitScoreCalculator.GetPearsonCorrelation(profile1, profile2) +
                    FitScoreCalculator.GetPearsonCorrelation(profile1, profile3));
        }        

        private double CalculateXicCorrelationOverTimeBetweenCharges(ChargeLcScanCluster cluster)
        {
            var maxCol = cluster.Members.Select(m => m.Col).Max();
            var minCol = cluster.Members.Select(m => m.Col).Min();

            var colLen = maxCol - minCol + 1;
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
            }
            
            var maxRow = cluster.Members.Select(m => m.Row).Max();
            var minRow = cluster.Members.Select(m => m.Row).Min();
            
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
            //var row1 = cluster.Members[0].Row;
            //var row2 = row1 + 1;
            //if (row2 == _intensityMap.Length) row2 = row1 - 1;

            //var profile1 = new double[maxCol - minCol + 1];
            //var profile2 = new double[maxCol - minCol + 1];
            //Array.Copy(_intensityMap[row1], minCol, profile1, 0, maxCol - minCol + 1);
            //Array.Copy(_intensityMap[row2], minCol, profile2, 0, maxCol - minCol + 1);

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
        
        public void BuildMatrix(double proteinMass)
        {
            const double epsillon = 1E-6;

            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            /*
            var isotopeIndexList = new List<int>();
            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                if (isoEnv.Envolope[isotopeIndex] > 0.0) isotopeIndexList.Add(isotopeIndex);

            _envelopePdf = new double[isotopeIndexList.Count];
            for (var i = 0; i < isotopeIndexList.Count; i++) _envelopePdf[i] = isoEnv.Envolope[isotopeIndexList[i]];

            var sum = _envelopePdf.Sum();
            for (var i = 0; i < isotopeIndexList.Count; i++) _envelopePdf[i] = _envelopePdf[i] / sum;
            */
            _envelopePdf        = isoEnv.Envolope;
            _intensityMap       = new double[_nCharges][];
            _correlationMap     = new double[_nCharges][];
            _xicMatrix          = new double[_nCharges][][];
            var observedCharges = new bool[_nCharges];

            
            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMap[i] = new double[_nScans];
                _correlationMap[i] = new double[_nScans];
                //_xicMatrix[i] = new double[_nScans][];
                //for (var j = 0; j < _nScans; j++) _xicMatrix[i][j] = new double[_envelopePdf.Length];
            }

            // Fill _intensityMap and _xicMatrix arrays
            //for (var i = 0; i < isotopeIndexList.Count; i++)
            for (var isotopeIndex = 0; isotopeIndex < _envelopePdf.Length; isotopeIndex++)
            {
                //var isotopeIndex = isotopeIndexList[i];
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    
                    // reading from file directly
                    //if (abundantIsotopeMz > _run.MaxMs1Mz || abundantIsotopeMz < _run.MinMs1Mz) continue;
                    //var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(abundantIsotopeMz - tolerance.GetToleranceAsTh(abundantIsotopeMz), abundantIsotopeMz + tolerance.GetToleranceAsTh(abundantIsotopeMz));
                    
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);
                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum) continue;
                    var xic = _cachedXic[curBinNum - _cachedMinBinNum];

                    if (!observedCharges[charge - _minCharge] && xic.Any(x => x > 0)) observedCharges[charge - _minCharge] = true;
                    //_intensityMap[charge - _minCharge] = _intensityMap[charge - _minCharge].Zip(xic, (a, b) => a + b).ToArray();

                    for (var j = 0; j < _nScans; j++)
                    {
                        _intensityMap[charge - _minCharge][j] += xic[j];
                        //_xicMatrix[charge - _minCharge][j][isotopeIndex] = xic[j];
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
                _xicMatrix[i] = new double[_nScans][];
                for (var j = 0; j < _nScans; j++)
                {
                    if (_intensityMap[i][j] < epsillon) continue;
                    _xicMatrix[i][j] = new double[_envelopePdf.Length];
                }

                var charge = i + _minCharge;
                for (var isotopeIndex = 0; isotopeIndex < _envelopePdf.Length; isotopeIndex++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);
                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum) continue;
                    var xic = _cachedXic[curBinNum - _cachedMinBinNum];

                    for (var j = 0; j < _nScans; j++)
                    {
                        if (_intensityMap[i][j] < epsillon) continue;

                        _xicMatrix[i][j][isotopeIndex] = xic[j];
                    }
                }
                
                for (var j = 0; j < _nScans; j++)
                {
                    if (_intensityMap[i][j] < epsillon) continue;
                    _correlationMap[i][j] = FitScoreCalculator.GetPearsonCorrelation(_xicMatrix[i][j], _envelopePdf);
                }
            }
        }

       
        private double GetIntensityThreshold(double topRatio)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            foreach (var t in _intensityMap)
                nonZeroIntensities.AddRange(t.Where(t1 => t1 > 0));

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[(int)Math.Floor(nonZeroIntensities.Count * (1-topRatio))];

            return intensityThreshold;
        }

        /*
        private List<ChargeScanCluster> FindClustersByExpansion(double scoreThreshold = 0)
        {
            var searchDirections = new ChargeScanCluster.ExpandDirection[4]
            {
                ChargeScanCluster.ExpandDirection.PrevScan, ChargeScanCluster.ExpandDirection.NextScan,
                ChargeScanCluster.ExpandDirection.ChargeUp, ChargeScanCluster.ExpandDirection.ChargeDown,
            };          
            
            var chargeObserved = new bool[_nCharges];
            var checkedOut = new bool[_nCharges, _nScans];

            for (var i = 0; i < _nCharges; i++)
            {
                if (_intensityMap[i].Sum() > 0) chargeObserved[i] = true;
            }

            var orderedSeedCells = GetSeedCellsByCorrelation(chargeObserved);
            var clusters = new List<ChargeScanCluster>();

            foreach (var seed in orderedSeedCells)
            {
                if (checkedOut[seed.Value.Row, seed.Value.Col]) continue;

                var newCluster = new ChargeScanCluster(_intensityMap, _data, seed.Value, seed.Key);
                var seedScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, newCluster.ObservedEnvelope);

                if (seedScore <= 0) continue;

                newCluster.SetScore(seedScore);

                for (var i = 0; i < searchDirections.Length;)
                {
                    var newEnvelop = newCluster.GetExpandedIsotopeEnvelop(searchDirections[i], chargeObserved);

                    if (newEnvelop == null)
                    {
                        i++;
                    }
                    else
                    {
                        var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, newEnvelop);
                        if (newScore >= newCluster.GetScore() * 0.9) newCluster.Expand(searchDirections[i], newEnvelop, newScore, ref checkedOut);
                        else i++;
                    }
                }
                
                if (newCluster.GetScore() > scoreThreshold) clusters.Add(newCluster);
            }

            return clusters.OrderByDescending(x => x.GetScore()).ToList();
        }*/
        

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
            const double corrLowerBound = 0.5;

            var checkedOut = new bool[_nCharges][];
            for (var i = 0; i < _nCharges; i++) checkedOut[i] = new bool[_nScans];
            
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seedCell in GetSeedCells(corrLowerBound))
            {
                if (checkedOut[seedCell.Row][seedCell.Col]) continue;

                var newCluster = new ChargeLcScanCluster(seedCell, _xicMatrix[seedCell.Row][seedCell.Col], _intensityMap[seedCell.Row][seedCell.Col], _correlationMap[seedCell.Row][seedCell.Col]);

                var neighbors = new Queue<ChargeLcScanCell>();
                neighbors.Enqueue(seedCell); // pick a seed
                checkedOut[seedCell.Row][seedCell.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();

                    for (var k = cell.Row - 2; k <= cell.Row + 2; k++)
                    {
                        if (k < _chargeIndexes.First() || k > _chargeIndexes.Last()) continue;

                        for (var l = cell.Col - 2; l <= cell.Col + 2; l++)
                        {
                            if (l < 0 || l >= _nScans || checkedOut[k][l] || _correlationMap[k][l] < corrLowerBound) continue;

                            //var tempEnvelope = _xicMatrix[k][l].Zip(newCluster.ObservedEnvelope, (a, b) => a + b).ToArray();
                            var tempEnvelope = new double[_envelopePdf.Length];
                            for (var t = 0; t < tempEnvelope.Length; t++) tempEnvelope[t] = newCluster.ObservedEnvelope[t] + _xicMatrix[k][l][t];

                            var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, tempEnvelope);

                            if (newScore >= newCluster.GetScore() || _correlationMap[k][l] > 0.85)
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

                if (newCluster.GetScore() > envelopCorrTh)
                {
                    var score2 = CalculateXicCorrelationOverTimeBetweenCharges(newCluster);
                    if (score2 > chargeCorrTh || newCluster.GetScore() > 0.9) clusters.Add(newCluster);

                    for (var i = newCluster.MinRow; i <= newCluster.MaxRow; i++)
                        for (var j = newCluster.MinCol; j <= newCluster.MaxCol; j++) checkedOut[i][j] = true;
                }
            }
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

    class ChargeLcScanCluster
    {
        private double _score;
        internal readonly List<ChargeLcScanCell> Members;
        internal double HighestIntensity { get; private set; }
        internal double[] ObservedEnvelope { get; private set; }
        internal double SeedCorrelation { get; private set; }

        internal int MaxCol { get; private set; }
        internal int MinCol { get; private set; }
        internal int MaxRow { get; private set; }
        internal int MinRow { get; private set; }

        internal ChargeLcScanCluster(ChargeLcScanCell seed, double[] seedEnvelope, double seedIntensity, double seedScore)
        {
            _score = seedScore;
            SeedCorrelation = seedScore;
            Members = new List<ChargeLcScanCell> {seed};
            ObservedEnvelope = seedEnvelope;
            HighestIntensity = seedIntensity;

            MaxCol = seed.Col;
            MinCol = seed.Col;
            MaxRow = seed.Row;
            MinRow = seed.Row;
        }

        internal ChargeScanRange GetChargeScanRange(int[] ms1ScanNums, int minCharge)
        {
            return new ChargeScanRange(MinRow + minCharge, MaxRow + minCharge, ms1ScanNums[MinCol], ms1ScanNums[MaxCol]);
        }

        internal double GetScore()
        {
            return _score;
        }

        internal void AddMember(ChargeLcScanCell member, double memberIntensity, double[] newObservedEnvelope, double newScore)
        {
            if (memberIntensity > HighestIntensity) HighestIntensity = memberIntensity;

            if (member.Col > MaxCol) MaxCol = member.Col;
            else if (member.Col < MinCol) MinCol = member.Col;
            if (member.Row > MaxRow) MaxRow = member.Row;
            else if (member.Row < MinRow) MinRow = member.Row;

            _score = newScore;
            ObservedEnvelope = newObservedEnvelope;
            Members.Add(member);
        }

        internal static int CompareByCorrelation(ChargeLcScanCluster x, ChargeLcScanCluster y)
        {
            return -x._score.CompareTo(y._score);
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
    /*
       internal class ChargeScanCluster
       {
           private double _score;
           internal readonly double HighestIntensity;
           internal double[] ObservedEnvelope { get; private set; }
           internal int MinChargeIndex { get; private set; }
           internal int MaxChargeIndex { get; private set; }
           internal int MinScanIndex { get; private set; }
           internal int MaxScanIndex { get; private set; }
           internal enum ExpandDirection {ChargeUp, ChargeDown, PrevScan, NextScan};

           private readonly double[,][] _data;
           private readonly double[][] _intensityMap;

           internal ChargeScanCluster(double[][] intensityMap, double[,][] data, ChargeLcScanCell seed, double seedIntensity)
           {
               _data = data;
               _intensityMap = intensityMap;
               _score = -1;
               MinChargeIndex = seed.Row;
               MaxChargeIndex = seed.Row;
               MinScanIndex = seed.Col;
               MaxScanIndex = seed.Col;
               HighestIntensity = seedIntensity;

               ObservedEnvelope = new double[_data.GetLength(1)];
               for (var i = 0; i < _data.GetLength(1); i++) ObservedEnvelope[i] = _data[seed.Row, i][seed.Col];                    
           }

           internal ChargeScanRange GetChargeScanRange(int[] ms1ScanNums, int minCharge)
           {
               return new ChargeScanRange(MinChargeIndex + minCharge, MaxChargeIndex + minCharge, ms1ScanNums[MinScanIndex], ms1ScanNums[MaxScanIndex]);
           }
        
           internal double GetScore()
           {
               return _score;
           }

           internal void SetScore(double score)
           {
               _score = score;
           }
           // if it's not expandable to the given direction, return null
           internal double[] GetExpandedIsotopeEnvelop(ExpandDirection direction, bool[] chargeObserved)
           {
               double[] newObservedEnvelope = null;

               if (direction == ExpandDirection.ChargeUp)
               {
                   if (MaxChargeIndex + 1 >= _data.GetLength(0) || chargeObserved[MaxChargeIndex + 1] == false) return null;

                   for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++)
                       if (_intensityMap[MaxChargeIndex + 1][scanIdx] < HighestIntensity*0.1) return null;

                   newObservedEnvelope = new double[_data.GetLength(1)];
                   Array.Copy(ObservedEnvelope, newObservedEnvelope, ObservedEnvelope.Length);

                   for (var i = 0; i < _data.GetLength(1); i++)
                   {
                       for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++)
                       {
                           ObservedEnvelope[i] += _data[MaxChargeIndex + 1, i][scanIdx];
                       }                    
                   }
               }
               else if (direction == ExpandDirection.ChargeDown)
               {
                   if (MinChargeIndex <= 0 || chargeObserved[MinChargeIndex - 1] == false) return null;

                   for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++)
                       if (_intensityMap[MinChargeIndex - 1][scanIdx] < HighestIntensity * 0.1) return null;

                   newObservedEnvelope = new double[_data.GetLength(1)];
                   Array.Copy(ObservedEnvelope, newObservedEnvelope, ObservedEnvelope.Length);
                
                   for (var i = 0; i < _data.GetLength(1); i++)
                   {
                       for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++)
                       {
                           ObservedEnvelope[i] += _data[MinChargeIndex - 1, i][scanIdx];
                       }
                   }
               }
               else if (direction == ExpandDirection.PrevScan)
               {
                   if (MinScanIndex <= 0) return null;

                   for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++)
                       if (_intensityMap[chargeIdx][MinScanIndex - 1] < HighestIntensity * 0.1) return null;

                   newObservedEnvelope = new double[_data.GetLength(1)];
                   Array.Copy(ObservedEnvelope, newObservedEnvelope, ObservedEnvelope.Length);
                   for (var i = 0; i < _data.GetLength(1); i++)
                   {
                       for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++)
                       {
                           ObservedEnvelope[i] += _data[chargeIdx, i][MinScanIndex - 1];
                       }
                   }
               }
               else if (direction == ExpandDirection.NextScan)
               {
                   if (MaxScanIndex + 1 >= _data[0,0].Length) return null;

                   for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++)
                       if (_intensityMap[chargeIdx][MaxScanIndex + 1] < HighestIntensity * 0.1) return null;

                   newObservedEnvelope = new double[_data.GetLength(1)];
                   Array.Copy(ObservedEnvelope, newObservedEnvelope, ObservedEnvelope.Length);
                
                   for (var i = 0; i < _data.GetLength(1); i++)
                   {
                       for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++)
                       {
                           ObservedEnvelope[i] += _data[chargeIdx, i][MaxScanIndex + 1];
                       }
                   }
               }

               return newObservedEnvelope;
           }

           internal void Expand(ExpandDirection directoin, double[] newObservedEnvelope, double newScore, ref bool[,] checkedOut)
           {
               switch (directoin)
               {
                   case ExpandDirection.ChargeUp:
                       MaxChargeIndex++;
                       for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++) checkedOut[MaxChargeIndex, scanIdx] = true;
                       break;
                   case ExpandDirection.ChargeDown:
                       MinChargeIndex--;
                       for (var scanIdx = MinScanIndex; scanIdx <= MaxScanIndex; scanIdx++) checkedOut[MinChargeIndex, scanIdx] = true;
                       break;
                   case ExpandDirection.NextScan:
                       MaxScanIndex++;
                       for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++) checkedOut[chargeIdx, MaxScanIndex] = true;
                       break;
                   case ExpandDirection.PrevScan:
                       MinScanIndex--;
                       for (var chargeIdx = MinChargeIndex; chargeIdx <= MaxChargeIndex; chargeIdx++) checkedOut[chargeIdx, MinScanIndex] = true;
                       break;
                   default:
                       break;
               }

               _score = newScore;
               ObservedEnvelope = newObservedEnvelope;
           }
       }
       */

}
