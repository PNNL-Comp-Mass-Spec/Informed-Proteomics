using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {

        public ChargeLcScanMatrix(LcMsRun run, int numBits = 26)
        {
            _minCharge = 2;
            _maxCharge = 50;
            _ms1ScanNums = run.GetMs1ScanVector();
            _numberOfScans = _ms1ScanNums.Length;
            _comparer = new MzComparerWithBinning(numBits);
            _nScans = _ms1ScanNums.Length;
            _nCharges = _maxCharge - _minCharge + 1;
            
            _run = run;
            _clusterMap = new Dictionary<int, IEnumerable<ChargeScanRange>>();
            _binNumToMs2ScanNumsMap = new Dictionary<int, IEnumerable<int>>();
            
            _cachedMinBinNum = _comparer.GetBinNumber(_run.MinMs1Mz);
            _cachedMaxBinNum = _comparer.GetBinNumber(_run.MaxMs1Mz);

            _cachedXic = new double[_cachedMaxBinNum - _cachedMinBinNum + 1][];

            for (var binNum = _cachedMinBinNum; binNum <= _cachedMaxBinNum; binNum++)
            {
                var mzStart = _comparer.GetMzStart(binNum);
                var mzEnd = _comparer.GetMzEnd(binNum);
                _cachedXic[binNum - _cachedMinBinNum] = _run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);
            }
        }

        public IEnumerable<ChargeScanRange> GetProbableChargeScanRegions(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);

            IEnumerable<ChargeScanRange> ranges;
            if (_clusterMap.TryGetValue(binNumber, out ranges)) return ranges;

            BuildMatrix(binNumber);
            var clusters = FindClusters(0.7);

            ranges = clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
            return _clusterMap[binNumber] = ranges;
        }

        public void PerformanceTest(double minMass, double maxMass)
        {
            var minBin = _comparer.GetBinNumber(minMass);
            var maxBin = _comparer.GetBinNumber(maxMass);

            for (var binNumber = minBin; binNumber <= maxBin; binNumber++)
            {
                BuildMatrix(binNumber);
                var clusters = FindClusters(0.7);
            }
        }
        
        public IEnumerable<int> GetMatchingMs2ScanNums(double monoIsotopicMass)
        {
            var binNumber = _comparer.GetBinNumber(monoIsotopicMass);
            IEnumerable<int> ms2Scans;
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

        public IEnumerable<ChargeScanRange> GetAllChargeScanRegions(double proteinMass, out double[] scores)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);
            var clusters = FindClusters(0.0);
            scores = clusters.Select(ac => ac.GetScore()).ToArray();
            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge));
        }

        private void BuildMatrix(int binNum)
        {
            double proteinMass = _comparer.GetMzAverage(binNum);
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            _envelopePdf        = new double[isoEnv.Envolope.Length];
            var zeroDummyXic    = new double[_nScans];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                _envelopePdf[isotopeIndex] = isoEnv.Envolope[isotopeIndex] / isoEnv.Envolope.Sum();

            _intensityMap = new double[_nCharges][];
            _data = new double[_nCharges, isoEnv.Envolope.Length][];

            for (var i = 0; i < _nCharges; i++)
            {
                _intensityMap[i] = new double[_numberOfScans];
            }

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
            {
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);

                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum)
                    {
                        _data[charge - _minCharge, isotopeIndex] = zeroDummyXic;
                    }
                    else
                    {
                        _intensityMap[charge - _minCharge] = _intensityMap[charge - _minCharge].Zip(_cachedXic[curBinNum - _cachedMinBinNum], (a, b) => a + b).ToArray();
                        _data[charge - _minCharge, isotopeIndex] = _cachedXic[curBinNum - _cachedMinBinNum];
                    }
                }
            }
        }
       

        private readonly LcMsRun _run;
        private readonly int _numberOfScans;
        private int _minCharge;
        private int _maxCharge;
        private int _nCharges;
        private readonly int _nScans;

        private double[,][] _data;
        private double[][] _intensityMap;

        private double[] _envelopePdf;

        private readonly int[] _ms1ScanNums;

        private readonly MzComparerWithBinning _comparer;

        // caching
        private readonly double[][] _cachedXic;
        private readonly int _cachedMinBinNum;
        private readonly int _cachedMaxBinNum;

        private readonly Dictionary<int, IEnumerable<ChargeScanRange>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;

        private double[] ExtractObservedProfile(int chargeIndex, int scanIndex, double[] initialValues = null)
        {
            var envelope = new double[_envelopePdf.Length];

            if (initialValues == null)
            {
                for (var i = 0; i < _envelopePdf.Length; i++) envelope[i] = _data[chargeIndex, i][scanIndex];
            }
            else
            {
                for (var i = 0; i < _envelopePdf.Length; i++) envelope[i] = _data[chargeIndex, i][scanIndex] + initialValues[i];
            }

            return envelope;
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
                intensityThreshold = nonZeroIntensities[nonZeroIntensities.Count / 2];

            return intensityThreshold;
        }
        
        private List<KeyValuePair<double, ChargeLcScanCell>> GetSeedCells(bool[] chargeObserved)
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
                    if (_intensityMap[i][j] > temp[j])
                    {
                        temp[j] = _intensityMap[i][j];
                        seedChargeIndex[j] = i;
                    }
                }
            }

            for (var j = 0; j < _nScans; j++)
            {
                if (temp[j] > 0)
                {
                    var pair = new KeyValuePair<double, ChargeLcScanCell>(temp[j], new ChargeLcScanCell(seedChargeIndex[j], j));
                    seedCells.Add(pair);
                }
            }

            return seedCells.OrderByDescending(cell => cell.Key).ToList();
        }
        private List<KeyValuePair<double, ChargeLcScanCell>> GetSeedCellsByIntensity(bool[] chargeObserved)
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

        public double CalcCorrelation(int minChargeIndex, int maxChargeIndex, int minScanIndex, int maxScanIndex)
        {

            var tempEnvelope = new double[_envelopePdf.Length];

            for (var i = minChargeIndex; i <= maxChargeIndex; i++)
            {
                for (var j = minScanIndex; j <= maxScanIndex; j++)
                {
                    tempEnvelope = ExtractObservedProfile(i, j, tempEnvelope);
                }
            }
            
            return FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, tempEnvelope);

        }

        private List<ChargeLcScanCluster> FindClusters(double scoreThreshold = 0)
        {
            var chargeObserved = new bool[_nCharges];
            
            for (var i = 0; i < _nCharges; i++)
            {
                if (_intensityMap[i].Sum() > 0) chargeObserved[i] = true;
                else chargeObserved[i] = false;
            }
            
            var orderedSeedCells = GetSeedCells(chargeObserved);
            
            var checkedOut = new byte[_nCharges, _nScans];
            foreach (var seed in orderedSeedCells) checkedOut[seed.Value.Row, seed.Value.Col] = 1;


            //1) Label start as 1
            //2) Get neighbor cells based on start point ->insert into a list L
            //3) Label neighbor as 2
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seed in orderedSeedCells)
            {
                if (checkedOut[seed.Value.Row, seed.Value.Col] == 2) continue;

                var seedEnvelope = ExtractObservedProfile(seed.Value.Row, seed.Value.Col);
                var newCluster = new ChargeLcScanCluster(seed.Value, seedEnvelope, seed.Key, _envelopePdf.Length);

                var neighbors = new Queue<ChargeLcScanCell>();
                neighbors.Enqueue(seed.Value); // pick a seed
                checkedOut[seed.Value.Row, seed.Value.Col] = 2;

                while (neighbors.Count > 0)
                {
                    var cell = neighbors.Dequeue();
                    //collect neighbor cells around the focused pixel
                    for (var k = cell.Row - 1; k < cell.Row + 2; k++)
                    {
                        if (k < 0 || k >= _nCharges || !chargeObserved[k]) continue;

                        for (var l = cell.Col - 1; l < cell.Col + 2; l++)
                        {
                            if (l < 0 || l >= _nScans || checkedOut[k, l] == 2 || _intensityMap[k][l] < newCluster.HighestIntensity * 0.1) continue;

                            var tempEnvelope = ExtractObservedProfile(k, l, newCluster.ObservedEnvelope);

                            var newScore = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, tempEnvelope);
                            var merge = false;

                            if (newScore >= newCluster.GetScore() || newScore > 0.7)
                            {
                                var newMember = new ChargeLcScanCell(k, l);
                                neighbors.Enqueue(newMember);
                                newCluster.AddMember(newMember, tempEnvelope, newScore);
                                checkedOut[k, l] = 2;
                            }
                        }
                    }
                }
                

                if (newCluster.GetScore() > scoreThreshold) clusters.Add(newCluster);
            }
            
            return clusters.OrderByDescending(x => x.GetScore()).ToList();
        }

        private double CalculateCorrelationXic(int chargeIndex1, int chargeIndex2, int minScanIndex, int maxScanIndex)
        {
            var scanLen = maxScanIndex - minScanIndex + 1;
            var xic1 = new double[scanLen];
            var xic2 = new double[scanLen];

            Array.Copy(_intensityMap[chargeIndex1], minScanIndex, xic1, 0, scanLen);
            Array.Copy(_intensityMap[chargeIndex2], minScanIndex, xic2, 0, scanLen);

            return FitScoreCalculator.GetPearsonCorrelation(xic1, xic2);
        }

        private double CalculateKlDivergence(double[] observedEnvelope)
        {
            double score = 0.0;
            //const double minIntensity = 1E-6;
            const double minIntensity = 1;
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

    class ChargeLcScanCluster
    {
        private double _score;
        internal readonly List<ChargeLcScanCell> Members;
        internal readonly double HighestIntensity;

        internal double[] ObservedEnvelope { get; private set; }

        internal ChargeLcScanCluster(ChargeLcScanCell seed, double[] seedEnvelope, double seedIntensity, int nEnvelope)
        {
            _score = -1;
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
