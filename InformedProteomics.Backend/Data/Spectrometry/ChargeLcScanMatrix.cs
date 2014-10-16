using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
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
        }

        public IEnumerable<ChargeScanRange> GetProbableChargeScanRegions(double monoIsotopicMass)
        {
            BuildMatrixWithoutBinning(monoIsotopicMass);
            var clusters = FindClusters(0.7);
            var ranges = clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
            return ranges;
            /*
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);

            IEnumerable<ChargeScanRange> ranges;
            if (_clusterMap.TryGetValue(binNumber, out ranges)) return ranges;

            BuildMatrix(binNumber);
            var clusters = FindClusters(0.7);

            ranges = clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
            return _clusterMap[binNumber] = ranges;
            */
        }
        
        public IEnumerable<int> GetMatchingMs2ScanNums(double monoIsotopicMass)
        {
            //IEnumerable<int> ms2Scans;
            //var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            //if (_binNumToMs2ScanNumsMap.TryGetValue(binNumber, out ms2Scans)) return ms2Scans;

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

            //_binNumToMs2ScanNumsMap.Add(binNumber, ms2ScanList);
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
            //var binNumber = _comparer.GetBinNumber(proteinMass);
            //BuildMatrix(binNumber);

            if (_cachedXic != null)
                BuildMatrix(proteinMass);
            else
                BuildMatrixWithoutBinning(proteinMass);

            var clusters = FindClusters(0.0);
            scores = new double[4][];
            for (var i = 0; i < scores.Length; i++) scores[i] = new double[clusters.Count];
            
            for(var i = 0; i < clusters.Count; i++)
            {
                //var tempEnv = new double[_envelopePdf.Length];
                //var maxCol = clusters[i].Members.Select(m => m.Col).Max();
                //var minCol = clusters[i].Members.Select(m => m.Col).Min();
                //var maxRow = clusters[i].Members.Select(m => m.Row).Max();
                //var minRow = clusters[i].Members.Select(m => m.Row).Min();
                //for (var y = minRow; y <= maxRow; y++) for (var x = minCol; x <= maxCol; x++) tempEnv = tempEnv.Zip(_xicData[y][x], (d1, d2) => d1 + d2).ToArray();
                //scores[0][i] = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, tempEnv);
                //scores[1][i] = CalculateKlDivergence(tempEnv);

                scores[0][i] = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, clusters[i].ObservedEnvelope);
                scores[1][i] = CalculateKlDivergence(clusters[i].ObservedEnvelope);
                scores[2][i] = CalculateXicCorrelationOverTime(clusters[i]);
                scores[3][i] = clusters[i].HighestIntensity;
            }
            
            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }

        /*
        public IEnumerable<ChargeScanRange> GetAllChargeScanRegionsNew(double proteinMass, out double[][] scores)
        {
            //var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(proteinMass);
            //BuildMatrix(binNumber);

            var clusters = FindClustersByExpansion(0.0);

            scores = new double[4][];
            scores[0] = new double[clusters.Count];
            scores[1] = new double[clusters.Count];
            scores[2] = new double[clusters.Count];
            scores[3] = new double[clusters.Count];

            for (var i = 0; i < clusters.Count; i++)
            {
                scores[0][i] = FitScoreCalculator.GetPearsonCorrelation(clusters[i].ObservedEnvelope, _envelopePdf);
                scores[1][i] = CalculateKlDivergence(clusters[i].ObservedEnvelope); 
                scores[2][i] = CalculateXicCorrelationOverTime(clusters[i]);
                scores[3][i] = clusters[i].HighestIntensity;
            }

            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }
        */

        

        private readonly LcMsRun _run;
        private int _minCharge;
        private int _maxCharge;
        private int _nCharges;
        private readonly int _nScans;

        private double[][] _intensityMap;

        private double[][][] _xicData;

        private double[] _envelopePdf;

        private readonly int[] _ms1ScanNums;

        private readonly MzComparerWithBinning _comparer;

        // caching
        private double[][] _cachedXic;
        private int _cachedMinBinNum;
        private int _cachedMaxBinNum;

        private readonly Dictionary<int, IEnumerable<ChargeScanRange>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;

        private double CalculateXicCorrelationOverTime(ChargeScanCluster cluster)
        {
            var row1 = (int)Math.Floor((cluster.MaxChargeIndex + cluster.MinChargeIndex) * 0.5);
            var row2 = row1 + 1;

            var profile1 = new double[cluster.MaxScanIndex - cluster.MinScanIndex + 1];
            var profile2 = new double[cluster.MaxScanIndex - cluster.MinScanIndex + 1];

            Array.Copy(_intensityMap[row1], cluster.MinScanIndex, profile1, 0, cluster.MaxScanIndex - cluster.MinScanIndex + 1);
            Array.Copy(_intensityMap[row2], cluster.MinScanIndex, profile2, 0, cluster.MaxScanIndex - cluster.MinScanIndex + 1);

            return FitScoreCalculator.GetPearsonCorrelation(profile1, profile2);
        }

        private double CalculateXicCorrelationOverTime(ChargeLcScanCluster cluster)
        {
            var maxCol = cluster.Members.Select(m => m.Col).Max();
            var minCol = cluster.Members.Select(m => m.Col).Min();
            var maxRow = cluster.Members.Select(m => m.Row).Max();
            var minRow = cluster.Members.Select(m => m.Row).Min();

            var row1 = (int)Math.Floor((maxRow + minRow) * 0.5);
            var row2 = row1 + 1;

            if (row2 == _intensityMap.Length)
            {
                row1--;
                row2--;
            }

            var profile1 = new double[maxCol - minCol + 1];
            var profile2 = new double[maxCol - minCol + 1];

            Array.Copy(_intensityMap[row1], minCol, profile1, 0, maxCol - minCol + 1);
            Array.Copy(_intensityMap[row2], minCol, profile2, 0, maxCol - minCol + 1);

            return FitScoreCalculator.GetPearsonCorrelation(profile1, profile2);
        }

        private void BuildMatrix(int binNum)
        {
            double proteinMass = _comparer.GetMzAverage(binNum);
            BuildMatrix(proteinMass);
        }

        public void BuildMatrixWithoutBinning(double proteinMass)
        {
            var tolerance = new Tolerance(10);
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            var isotopeIndexList = new List<int>();
            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                if (isoEnv.Envolope[isotopeIndex] > 0.0) isotopeIndexList.Add(isotopeIndex);

            _envelopePdf = new double[isotopeIndexList.Count];
            for (var i = 0; i < isotopeIndexList.Count; i++) _envelopePdf[i] = isoEnv.Envolope[isotopeIndexList[i]];

            var sum = _envelopePdf.Sum();
            for (var i = 0; i < isotopeIndexList.Count; i++) _envelopePdf[i] = _envelopePdf[i] / sum;

            _intensityMap = new double[_nCharges][];

            _xicData = new double[_nCharges][][];
            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMap[i] = new double[_nScans];
                _xicData[i] = new double[_nScans][];
                for (var j = 0; j < _nScans; j++)
                {
                    _xicData[i][j] = new double[_envelopePdf.Length];
                }
            }

            for (var i = 0; i < isotopeIndexList.Count; i++)
            {
                var isotopeIndex = isotopeIndexList[i];
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);

                    if (abundantIsotopeMz > _run.MaxMs1Mz || abundantIsotopeMz < _run.MinMs1Mz) continue;

                    var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(abundantIsotopeMz - tolerance.GetToleranceAsTh(abundantIsotopeMz), abundantIsotopeMz + tolerance.GetToleranceAsTh(abundantIsotopeMz));

                    _intensityMap[charge - _minCharge] = _intensityMap[charge - _minCharge].Zip(xic, (a, b) => a + b).ToArray();

                    for (var j = 0; j < _nScans; j++) _xicData[charge - _minCharge][j][i] = xic[j];
                }
            }
        }

        private void BuildMatrix(double proteinMass)
        {
            //double proteinMass = _comparer.GetMzAverage(binNum);
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);
            var zeroDummyXic = new double[_nScans];

            var isotopeIndexList = new List<int>();
            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                if (isoEnv.Envolope[isotopeIndex] > 0.0) isotopeIndexList.Add(isotopeIndex);

            _envelopePdf = new double[isotopeIndexList.Count];
            for (var i = 0; i < isotopeIndexList.Count; i++)
            {
                var isotopeIndex = isotopeIndexList[i];
                _envelopePdf[i] = isoEnv.Envolope[isotopeIndex];
            }

            var sum = _envelopePdf.Sum();
            for (var i = 0; i < isotopeIndexList.Count; i++) _envelopePdf[i] = _envelopePdf[i] / sum;

            _intensityMap = new double[_nCharges][];
            _xicData = new double[_nCharges][][];
            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMap[i] = new double[_nScans];
                _xicData[i] = new double[_nScans][];
                for (var j = 0; j < _nScans; j++)
                {
                    _xicData[i][j] = new double[_envelopePdf.Length];
                }
            }

            //for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
            for (var i = 0; i < isotopeIndexList.Count; i++)
            {
                var isotopeIndex = isotopeIndexList[i];

                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);

                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum) continue;

                    _intensityMap[charge - _minCharge] = _intensityMap[charge - _minCharge].Zip(_cachedXic[curBinNum - _cachedMinBinNum], (a, b) => a + b).ToArray();

                    for (var j = 0; j < _nScans; j++) _xicData[charge - _minCharge][j][i] = _cachedXic[curBinNum - _cachedMinBinNum][j];
                }
            }
        }       

        private double GetIntensityThreshold(double[,] matrix)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            for (var i = 0; i < matrix.GetLength(0); i++)
                for (var j = 0; j < matrix.GetLength(1); j++)
                    if (matrix[i, j] > 0) nonZeroIntensities.Add(matrix[i, j]);

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[nonZeroIntensities.Count / 2];

            return intensityThreshold;
        }

        private double GetIntensityThreshold(IEnumerable<double[]> matrix)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            foreach (var t in matrix)
                nonZeroIntensities.AddRange(t.Where(t1 => t1 > 0));

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[nonZeroIntensities.Count/2];

            return intensityThreshold;
        }

        private IEnumerable<KeyValuePair<double, ChargeLcScanCell>> GetSeedCellsByCorrelation(bool[] chargeObserved)
        {
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();
            var intensityThreshold = GetIntensityThreshold(_intensityMap);

            for (var i = 0; i < _nCharges; i++)
            {
                if (!chargeObserved[i]) continue;

                for (var j = 0; j < _nScans; j++)
                {
                    if (!(_intensityMap[i][j] > intensityThreshold)) continue; 
                    var corr = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, _xicData[i][j]);

                    if (corr <= 10E-3) continue; 
                    var pair = new KeyValuePair<double, ChargeLcScanCell>(corr, new ChargeLcScanCell(i, j));
                    seedCells.Add(pair);
                }
            }

            return seedCells.OrderByDescending(cell => cell.Key);
        }

        private IEnumerable<KeyValuePair<double, ChargeLcScanCell>> GetSeedCellsForEveryScans(bool[] chargeObserved)
        {
            // Pick the max intensity charge for each scan time
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();
            var seedChargeIndex = new int[_nScans];
            var temp = new double[_nScans];
            
            for (var i = 0; i < _nCharges; i++ )
            {
                if (!chargeObserved[i]) continue;

                for (var j = 0; j < _nScans; j++)
                {
                    if (!(_intensityMap[i][j] > temp[j])) continue;
                    temp[j] = _intensityMap[i][j];
                    seedChargeIndex[j] = i;
                }
            }

            var intensityThreshold = GetIntensityThreshold(_intensityMap);
            for (var j = 0; j < _nScans; j++)
            {
                if (temp[j] <= intensityThreshold) continue;
                var pair = new KeyValuePair<double, ChargeLcScanCell>(temp[j], new ChargeLcScanCell(seedChargeIndex[j], j));
                seedCells.Add(pair);
            }

            return seedCells.OrderByDescending(cell => cell.Key).ToList();
        }

        private IEnumerable<KeyValuePair<double, ChargeLcScanCell>> GetSeedCellsByIntensity(bool[] chargeObserved)
        {
            // Pick cells with intensities over a certain threshold 
            var seedCells = new List<KeyValuePair<double, ChargeLcScanCell>>();
            var intensityThreshold = GetIntensityThreshold(_intensityMap);
            for (var i = 0; i < _nCharges; i++)
            {
                if (!chargeObserved[i]) continue;
                for (var j = 0; j < _nScans; j++)
                {
                    if (_intensityMap[i][j] > intensityThreshold)
                    {
                        var pair = new KeyValuePair<double, ChargeLcScanCell>(_intensityMap[i][j], new ChargeLcScanCell(i, j));
                        seedCells.Add(pair);
                    }
                }
            }
            return seedCells.OrderByDescending(cell => cell.Key).ToList();
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

        private List<ChargeLcScanCluster> FindClusters(double scoreThreshold = 0)
        {
            var chargeObserved  = new bool[_nCharges];
            var checkedOut      = new bool[_nCharges][];

            for (var i = 0; i < _nCharges; i++)
            {
                if (_intensityMap[i].Sum() > 0) chargeObserved[i] = true;
                checkedOut[i] = new bool[_nScans];
            }

            var orderedSeedCells = GetSeedCellsForEveryScans(chargeObserved);
            
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seed in orderedSeedCells)
            {
                if (checkedOut[seed.Value.Row][seed.Value.Col]) continue;

                var seedCell = seed.Value;
                var seedScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf,_xicData[seedCell.Row][seedCell.Col]);

                if (seedScore < scoreThreshold) continue;

                var newCluster = new ChargeLcScanCluster(seedCell, _xicData[seedCell.Row][seedCell.Col], _intensityMap[seedCell.Row][seedCell.Col], seedScore);

                var neighbors = new Queue<ChargeLcScanCell>();
                neighbors.Enqueue(seed.Value); // pick a seed
                checkedOut[seed.Value.Row][seed.Value.Col] = true;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();

                    //collect neighbor cells around the focused pixel
                    for (var k = cell.Row - 1; k < cell.Row + 2; k++)
                    {
                        if (k < 0 || k >= _nCharges || !chargeObserved[k]) continue;

                        for (var l = cell.Col - 1; l < cell.Col + 2; l++)
                        {
                            if (l < 0 || l >= _nScans || checkedOut[k][l] || _intensityMap[k][l] < newCluster.HighestIntensity * 0.1) continue;

                            var tempEnvelope = _xicData[k][l].Zip(newCluster.ObservedEnvelope, (a, b) => a + b).ToArray();
                            var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, tempEnvelope);

                            if (!(newScore >= newCluster.GetScore())) continue;

                            var newMember = new ChargeLcScanCell(k, l);
                            neighbors.Enqueue(newMember);
                            newCluster.AddMember(newMember, tempEnvelope, newScore);
                            checkedOut[k][l] = true;
                        }
                    }
                }

                
                if (newCluster.GetScore() > scoreThreshold && CalculateKlDivergence(newCluster.ObservedEnvelope) < 0.3) clusters.Add(newCluster);
            }
            
            return clusters.OrderByDescending(x => x.GetScore()).ToList();
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
                var xic = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);

                _cachedXic[binNum - _cachedMinBinNum] = xic;
            }
        }
            
        public void PrepareSmoothFullPrecursorXicVector()
        {
            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);
            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];

            var ms1ToIndex = new int[_ms1ScanNums[_ms1ScanNums.Length - 1] + 1];
            for (var i = 0; i < _ms1ScanNums.Length; i++)
            {
                ms1ToIndex[_ms1ScanNums[i]] = i;
            }

            var binNumToAvgMass = new double[_cachedMaxBinNum - _cachedMinBinNum + 1];
            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                _cachedXic[binNum - _cachedMinBinNum] = new double[_nScans];
                binNumToAvgMass[binNum - _cachedMinBinNum] = _comparer.GetMzAverage(binNum);
            }

            var xic = ((PbfLcMsRun)_run).GetAllPrecursorExtractedIonChromatogram();
            var tolerance = _comparer.GetTolerance();
            
            const int binNumWindow = 2;
            foreach (var xicPoint in xic)
            {
                var xicBinNum = _comparer.GetBinNumber(xicPoint.Mz);
                var bandWidth = tolerance.GetToleranceAsTh(xicPoint.Mz) * 0.5;
                //var bandWidth = tolerance.GetToleranceAsTh(xicPoint.Mz);
                
                for (var binNum = Math.Max(xicBinNum - binNumWindow, _cachedMinBinNum); binNum <= Math.Min(xicBinNum + binNumWindow, _cachedMaxBinNum); binNum++)
                {
                    var u = (xicPoint.Mz - binNumToAvgMass[binNum - _cachedMinBinNum]) / bandWidth;
                    var smoothedIntensity = xicPoint.Intensity*(Math.Exp(-Math.Pow(u, 2)/2));
                    _cachedXic[binNum - _cachedMinBinNum][ms1ToIndex[xicPoint.ScanNum]] += smoothedIntensity;
                }
            }
        }
        * /
        
        /*private double CalculateCorrelationXic(int chargeIndex1, int chargeIndex2, int minScanIndex, int maxScanIndex)
        {
            var scanLen = maxScanIndex - minScanIndex + 1;
            var xic1 = new double[scanLen];
            var xic2 = new double[scanLen];

            Array.Copy(_intensityMap[chargeIndex1], minScanIndex, xic1, 0, scanLen);
            Array.Copy(_intensityMap[chargeIndex2], minScanIndex, xic2, 0, scanLen);

            return FitScoreCalculator.GetPearsonCorrelation(xic1, xic2);
        }*/

        private double CalculateKlDivergence(double[] observedEnvelope)
        {
            double score = 0.0;
            const double minIntensity = 1E-9;
            //const double minIntensity = 1;
            var normalizedEnv = new double[observedEnvelope.Length];

            var sum = observedEnvelope.Sum();
            bool foundZero = false;
            for (var k = 0; k < _envelopePdf.Length; k++)
            {
                normalizedEnv[k] = observedEnvelope[k] / sum;
                if (!(normalizedEnv[k] < minIntensity)) continue;
                normalizedEnv[k] = minIntensity;
                foundZero = true;
            }

            if (!foundZero) return SimpleMath.GetKLDivergence(_envelopePdf, normalizedEnv);
            sum = normalizedEnv.Sum();
            for (var k = 0; k < _envelopePdf.Length; k++) normalizedEnv[k] = normalizedEnv[k] / sum;

            score = SimpleMath.GetKLDivergence(_envelopePdf, normalizedEnv);
            
            return score;
        }
   }

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


    class ChargeLcScanCluster
    {
        private double _score;
        internal readonly List<ChargeLcScanCell> Members;
        internal readonly double HighestIntensity;

        internal double[] ObservedEnvelope { get; private set; }

        internal ChargeLcScanCluster(ChargeLcScanCell seed, double[] seedEnvelope, double seedIntensity, double seedScore)
        {
            _score = seedScore;
            Members = new List<ChargeLcScanCell> {seed};
            ObservedEnvelope = seedEnvelope;
            HighestIntensity = seedIntensity;
        }

        internal ChargeScanRange GetChargeScanRange(int[] ms1ScanNums, int minCharge)
        {
            var maxCol = Members.Select(m => m.Col).Max();
            var minCol = Members.Select(m => m.Col).Min();
            var maxRow = Members.Select(m => m.Row).Max();
            var minRow = Members.Select(m => m.Row).Min();

            return new ChargeScanRange(minRow + minCharge, maxRow + minCharge, ms1ScanNums[minCol], ms1ScanNums[maxCol]);
        }

        internal double GetScore()
        {
            return _score;
        }

        internal void AddMember(ChargeLcScanCell member, double[] newObservedEnvelope, double newScore)
        {
            Members.Add(member);
            _score = newScore;
            ObservedEnvelope = newObservedEnvelope;
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
