using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    class RankProbability
    {
        public Dictionary<IonType, int> IonFrequencies { get; private set; }
        public int RankCount { get; set; }

        public Dictionary<IonType, IonProbability> IonProbabilities
        {
            get
            {
                var probDict = IonFrequencies.Keys.ToDictionary(ionType => ionType, ionType => new IonProbability(ionType.Name, IonFrequencies[ionType], RankCount));
                return probDict;
            }
        }

        public RankProbability(IEnumerable<IonType> ionTypes)
        {
            IonFrequencies = new Dictionary<IonType, int>();
            RankCount = 0;

            foreach (var ionType in ionTypes)
            {
                IonFrequencies.Add(ionType, 0);
            }
        }
    }

    public class RankTable
    {
        private readonly IonType[] _ionTypes;
        private readonly List<RankProbability> _rankTable;

        public int TotalRanks { get { return _rankTable.Count;  } }

        public RankTable(IonType[] ionTypes)
        {
            _ionTypes = ionTypes;
            _rankTable = new List<RankProbability>();
        }

        public List<Dictionary<IonType, IonProbability>> IonProbabilities
        {
            get
            {
                return _rankTable.Select(rankProbability => rankProbability.IonProbabilities).ToList();
            }
        }

        public IonProbability[,] RankProbabilities
        {
            get
            {
                var ranks = new IonProbability[_rankTable.Count, _ionTypes.Length];
                for (int i = 0; i < _rankTable.Count; i++)
                {
                    for (int j = 0; j < _ionTypes.Length; j++)
                    {
                        var ionType = _ionTypes[j];
                        ranks[i, j] = new IonProbability(ionType.Name, _rankTable[i].IonFrequencies[ionType], _rankTable[i].RankCount);
                    }
                }
                return ranks;
            }
        }

        private void ExtendRankTable(int size)
        {
            while (_rankTable.Count < size)
            {
                _rankTable.Add(new RankProbability(_ionTypes));
            }
        }

        public void RankMatches(IEnumerable<SpectrumMatch> matches, Tolerance tolerance)
        {
            foreach (var match in matches)
            {
                var ranks = new RankedPeaks(match.Spectrum);
                foreach (var ionType in _ionTypes)
                {
                    var ions = match.GetCleavageIons(ionType);
                    foreach (var ion in ions)
                    {
                        ranks.RankIon(ionType, ion, tolerance);
                    }
                }
                AddRanks(ranks);
            }
        }

        public void AddRanks(RankedPeaks rankedPeaks)
        {
            var ranks = rankedPeaks.Ranks;
            ExtendRankTable(ranks.Count);

            for (int i = 0; i < ranks.Count; i++)
            {
                _rankTable[i].RankCount++;
                if (ranks[i].Iontype != null)
                    _rankTable[i].IonFrequencies[ranks[i].Iontype]++;
            }
        }
    }
}
