using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeLcScanMatrix : ILcMsMap, ISequenceFilter
    {
        int NRows { get { return _maxCharge - _minCharge + 1; } }
        int NColumns { get { return _numberOfScans; } }
        
        public enum ScoringMeasure
        {
            PearsonCorrelation,
            KlDivergence,
            Hybrid
        };

        public ChargeLcScanMatrix(LcMsRun run, int numBits = 26)
        {
            _minCharge = 2;
            _maxCharge = 60;
            _ms1ScanNums = run.GetMs1ScanVector();
            _numberOfScans = _ms1ScanNums.Length;
            _comparer = new MzComparerWithBinning(numBits);
            
            _run = run;
            _clusterMap = new Dictionary<int, IEnumerable<ChargeScanRange>>();
            _binNumToMs2ScanNumsMap = new Dictionary<int, IEnumerable<int>>();

            _scoringType = ScoringMeasure.PearsonCorrelation;
            _scoreThreshold = 0.7;
            
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
            var clusters = GenerateClusters(_sumOverIsotope);
            ranges = clusters.Where(cluster => cluster.GetScore() > _scoreThreshold)
                        .Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge)).ToList();
            _clusterMap[binNumber] = ranges;

            return ranges;
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
 
        public void SetScoringType(ScoringMeasure scoringType)
        {
            _scoringType = scoringType;

            if (_scoringType == ScoringMeasure.PearsonCorrelation) _scoreThreshold = 0.7;
            else if (_scoringType == ScoringMeasure.KlDivergence) _scoreThreshold = 0.5;
        }

        public void SetChargeRange(int minCharge, int maxCharge)
        {
            _minCharge = minCharge;
            _maxCharge = maxCharge;
        }

        public IEnumerable<ChargeScanRange> GetAllChargeScanRegions(double proteinMass, out double[] scores)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);
            var clusters = GenerateClusters(_sumOverIsotope);
            scores = clusters.Select(ac => ac.GetScore()).ToArray();
            return clusters.Select(cluster => cluster.GetChargeScanRange(_ms1ScanNums, _minCharge)).ToList();
        }


        private void BuildMatrix(int binNum)
        {
            double proteinMass = _comparer.GetMzAverage(binNum);
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            _envelopePdf = new double[isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
                _envelopePdf[isotopeIndex] = isoEnv.Envolope[isotopeIndex]/isoEnv.Envolope.Sum();
            
            _sumOverIsotope = new double[NRows, NColumns];
            _data = new double[NRows, NColumns, isoEnv.Envolope.Length];

            for (var isotopeIndex = 0; isotopeIndex < isoEnv.Envolope.Length; isotopeIndex++)
            {
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var abundantIsotopeMz = Ion.GetIsotopeMz(proteinMass, charge, isotopeIndex);
                    var curBinNum = _comparer.GetBinNumber(abundantIsotopeMz);

                    if (curBinNum < _cachedMinBinNum || curBinNum > _cachedMaxBinNum) continue;

                    for (var scanIndex = 0; scanIndex < _numberOfScans; scanIndex++)
                    {
                        _sumOverIsotope[charge - _minCharge, scanIndex] += _cachedXic[curBinNum - _cachedMinBinNum][scanIndex];
                        _data[charge - _minCharge, scanIndex, isotopeIndex] = _cachedXic[curBinNum - _cachedMinBinNum][scanIndex];
                    }
                }
            }
        }

        public string PrintMatrix(double proteinMass)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);

            var sb = new StringBuilder();

            for (var i = 0; i < NRows; i++)
            {
                for (var j = 0; j < NColumns; j++)
                {
                    sb.Append(_sumOverIsotope[i, j]);

                    if (j != NColumns - 1) sb.Append("\t");
                }
                sb.Append("\n");
            }

            return sb.ToString();
        }

        public double[, ,] GetRawMatrix(double proteinMass)
        {
            var binNumber = _comparer.GetBinNumber(proteinMass);
            BuildMatrix(binNumber);

            return _data; 
        }

        private readonly LcMsRun _run;
        private readonly int _numberOfScans;
        private int _minCharge;
        private int _maxCharge;
        private double[, ,] _data;
        private double[,] _sumOverIsotope;
        private double[] _envelopePdf;

        private double _scoreThreshold = 0.5;
        private readonly int[] _ms1ScanNums;

        private readonly MzComparerWithBinning _comparer;

        // caching
        private readonly double[][] _cachedXic;
        private readonly int _cachedMinBinNum;
        private readonly int _cachedMaxBinNum;

        private readonly Dictionary<int, IEnumerable<ChargeScanRange>> _clusterMap; // BinNum -> IEnumerable<ChargeScanRange>
        private readonly Dictionary<int, IEnumerable<int>> _binNumToMs2ScanNumsMap;

        private ScoringMeasure _scoringType;

        private double GetIntensityThreshold(double[,] matrix)
        {
            double intensityThreshold = 0;
            var nonZeroIntensities = new List<double>();

            for (var i = 0; i < NRows; i++)
                for (var j = 0; j < NColumns; j++)
                    if (matrix[i, j] > 0) nonZeroIntensities.Add(matrix[i, j]);

            nonZeroIntensities.Sort();

            if (nonZeroIntensities.Count > 2)
                intensityThreshold = nonZeroIntensities[nonZeroIntensities.Count / 2];

            return intensityThreshold;
        }

        private double ComputeSingleCellScore(ChargeLcScanCell newCell)
        {
            var tempEnvelope = new double[_envelopePdf.Length];

            for (var k = 0; k < _envelopePdf.Length; k++)
                tempEnvelope[k] = _data[newCell.Row, newCell.Col, k];

            return ComputeScore(tempEnvelope, _scoringType);            
        }

        private double ComputeClusterScoreWithCell(ChargeLcScanCluster cluster, ChargeLcScanCell newCell)
        {
            var tempEnvelope = new double[_envelopePdf.Length];
            for (var k = 0; k < _envelopePdf.Length; k++)
                tempEnvelope[k] = cluster.ObservedEnvelope[k] + _data[newCell.Row, newCell.Col, k];

            return ComputeScore(tempEnvelope, _scoringType);
        }
        
        private double ComputeScore(double[] observedEnvelope, ScoringMeasure scoringType)
        {
            double score = 0.0;

            if (scoringType == ScoringMeasure.PearsonCorrelation)
            {
                score = FitScoreCalculator.GetPearsonCorrelation(_envelopePdf, observedEnvelope);
            }
            else if (scoringType == ScoringMeasure.KlDivergence)
            {
                const double minIntensity = 1E-6;
                var normalizedEnv = new double[observedEnvelope.Length];

                var sum = observedEnvelope.Sum();
                bool foundZero = false;
                for (var k = 0; k < _envelopePdf.Length; k++)
                {
                    normalizedEnv[k] = observedEnvelope[k] / sum;
                    if (normalizedEnv[k] < minIntensity)
                    {
                        normalizedEnv[k] = minIntensity;
                        foundZero = true;
                    }
                }

                if (!foundZero) return SimpleMath.GetKLDivergence(_envelopePdf, normalizedEnv);

                sum = normalizedEnv.Sum();
                for (var k = 0; k < _envelopePdf.Length; k++)
                    normalizedEnv[k] = normalizedEnv[k] / sum;

                score = SimpleMath.GetKLDivergence(_envelopePdf, normalizedEnv);

            }
            else if (scoringType == ScoringMeasure.Hybrid)
            {
                
            }

            return score;
        }
        
        private IEnumerable<ChargeLcScanCluster> GenerateClusters(double[,] matrix)
        {
            var intensityThreshold = GetIntensityThreshold(matrix);
            var tempMap = new byte[NRows, NColumns];
            var seedCells = new List<ChargeLcScanCell>();

            
            for(var i = 0; i < NRows; i++)
                for (var j = 0; j < NColumns; j++)
                    if (matrix[i, j] > intensityThreshold)
                    {
                        tempMap[i, j] = 1;
                        seedCells.Add(new ChargeLcScanCell(i, j, matrix[i, j]));
                    }

            var orderedSeedCells = seedCells.OrderByDescending(cell => cell.Intensity);

            //1) Label start as 1
	        //2) Get neighbor cells based on start point ->insert into a list L
	        //3) Label neighbor as 2
            var clusters = new List<ChargeLcScanCluster>();

            foreach (var seed in orderedSeedCells)
            {
                if (tempMap[seed.Row, seed.Col] != 1) continue;

                var newCluster = new ChargeLcScanCluster(_envelopePdf.Length);
                var neighbors = new Queue<ChargeLcScanCell>();

                neighbors.Enqueue(seed); // pick a seed
                
                newCluster.Add(seed, _data, -1.0);
                tempMap[seed.Row, seed.Col] = 2;  

                do
                {
                    var cell = neighbors.Dequeue();
                    //collect neighbor cells around the focused pixel
                    for (var k = cell.Row - 1; k < cell.Row + 2; k++)
                    {
                        for (var l = cell.Col - 1; l < cell.Col + 2; l++)
                        {
                            if (k < 0 || l < 0 || k >= NRows || l >= NColumns || tempMap[k, l] == 2 || matrix[k, l] < newCluster.HighestIntensity * 0.1) continue;

                            var tempEnvelope = new double[_envelopePdf.Length];
                            for (var i = 0; i < _envelopePdf.Length; i++)
                                tempEnvelope[i] = newCluster.ObservedEnvelope[i] + _data[k, l, i];

                            var newScore = ComputeScore(tempEnvelope, _scoringType);
                            var merge = false;

                            if (_scoringType == ScoringMeasure.PearsonCorrelation && newScore >= newCluster.GetScore() || newScore > 0.7)
                            {
                                merge = true;
                            }
                            else if (_scoringType == ScoringMeasure.KlDivergence)
                            {
                                var pearsonCorr = ComputeScore(tempEnvelope, ScoringMeasure.PearsonCorrelation);
                                if (newScore <= newCluster.GetScore() || pearsonCorr > 0.7) merge = true;
                            }

                            if (!merge) continue;
                            
                            var newMember = new ChargeLcScanCell(k, l, matrix[k, l]);
                            neighbors.Enqueue(newMember);
                            newCluster.Add(newMember, _data, newScore);
                            tempMap[k, l] = 2;
                        }
                    }
                }
                while (neighbors.Count > 0);

                if (newCluster.Members.Count == 1)
                {
                    tempMap[newCluster.Members[0].Row, newCluster.Members[0].Col] = 1;
                }
                else clusters.Add(newCluster);                 
            }

            if (_scoringType == ScoringMeasure.PearsonCorrelation)
                return clusters.OrderByDescending(x => x.GetScore()).ToList();
            if (_scoringType == ScoringMeasure.KlDivergence)
                return clusters.OrderBy(x => x.GetScore()).ToList();

            return null;
        }
   }

    class ChargeLcScanCluster
    {
        private double _score;

        internal readonly List<ChargeLcScanCell> Members;
        internal double[] ObservedEnvelope { get; private set; }

        internal ChargeLcScanCluster(int nEnvelope)
        {
            _score = -1;
            Members = new List<ChargeLcScanCell>();
            ObservedEnvelope = new double[nEnvelope];
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

        internal void Add(ChargeLcScanCell cell, double[,,] data, double newScore)
        {
            Members.Add(cell);
            _score = newScore;
            for (var k = 0; k < ObservedEnvelope.Length; k++)
                ObservedEnvelope[k] += data[cell.Row, cell.Col, k];
        }

        internal double HighestIntensity
        {
            get { return Members.Count < 1 ? 0.0 : Members[0].Intensity; }
        }
    }

    class ChargeLcScanCell
    {
        internal int Row;
        internal int Col;
        internal double Intensity;

        internal ChargeLcScanCell(int row, int col, double intensity)
        {
            Row = row;
            Col = col;
            Intensity = intensity;
        }
    }
}
