using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class RankTable: I2DProbabilityTable<IonType>
    {
        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly List<RankProbability> _rankTable;

        public RankProbability NotFound { get; set; }
        public int MaxRanks { get; private set; }
        public int TotalRanks { get { return _rankTable.Count;  } }

        public RankTable(IonType[] ionTypes, Tolerance tolerance, int maxRanks)
        {
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            MaxRanks = maxRanks;
            _rankTable = new List<RankProbability> {Capacity = MaxRanks};
            for (int i = 0; i < MaxRanks; i++)
            {
                _rankTable.Add(new RankProbability(_ionTypes));
            }
            NotFound = new RankProbability(_ionTypes);
        }

        public RankTable(List<RankProbability> rankTable, IonType[] ionTypes, Tolerance tolerance, int maxRanks)
        {
            _rankTable = rankTable;
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            MaxRanks = maxRanks;
        }

        public List<Dictionary<IonType, Probability<string>>> IonProbabilities
        {
            get
            {
                var rankProbabilities = _rankTable.Select(rankProbability => rankProbability.IonProbabilities).ToList();
                rankProbabilities.Add(NotFound.IonProbabilities);
                return rankProbabilities;
            }
        }

        public void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                var ranks = new RankedPeaks(match.Spectrum);
                foreach (var ionType in _ionTypes)
                {
                    var ions = match.GetCleavageIons(ionType);
                    foreach (var ion in ions)
                    {
                        ranks.RankIon(ionType, ion, _tolerance);
                    }
                }
                AddRanks(ranks);
            }
        }

        public void Smooth(int[] smoothingRanks, int[] smoothingWindowSize)
        {
            if (smoothingRanks.Length != smoothingWindowSize.Length)
            {
                throw new ArgumentException("Unequal length smoothingRanks and smoothingWindowSize.");
            }

            foreach (var ionType in _ionTypes)
            {
                int index = 0;
                int startRank = 0;
                int window = smoothingWindowSize[index];
                var total = 0.0;
                var count = 0;

                for (int currRank = startRank; currRank < MaxRanks; currRank += window+1)
                {
                    var endRank = Math.Min(currRank + window + 1, MaxRanks);
                    for (int i = currRank; i < endRank; i++)
                    {
                        total += _rankTable[i].IonFrequencies[ionType];
                        count++;
                    }
                    var average = total/count;
                    if (!average.Equals(0))
                    {
                        for (int i = startRank; i < endRank; i++)
                        {
                            _rankTable[i].IonFrequencies[ionType] = average;
                        }
                        count = 0;
                        total = 0;
                        startRank = currRank;
                    }
                    if (currRank >= smoothingRanks[index])
                    {
                        index++;
                        window = smoothingWindowSize[index];
                    }
                    if ((currRank + 2 * window) > MaxRanks)
                    {
                        window *= 2;
                    }
                }
            }
        }

        public void AddRanks(RankedPeaks rankedPeaks)
        {
            var ranks = rankedPeaks.Ranks;

            for (int i = 0; (i < ranks.Count && i < MaxRanks-1); i++)
            {
                _rankTable[i].RankCount++;
                if (ranks[i].Iontypes.Count > 0)
                {
                    foreach (var ionType in ranks[i].Iontypes)
                    {
                        _rankTable[i].IonFrequencies[ionType]++;
                    }
                }
            }
            if (ranks.Count > MaxRanks - 1)
            {
                for (int i = MaxRanks - 1; i < ranks.Count; i++)
                {
                    _rankTable[MaxRanks - 1].RankCount++;
                    if (ranks[i].Iontypes.Count > 0)
                    {
                        foreach (var ionType in ranks[i].Iontypes)
                        {
                            _rankTable[MaxRanks - 1].IonFrequencies[ionType]++;
                        }
                    }
                }
            }
            NotFound.AddIons(rankedPeaks.NotFound);
            NotFound.RankCount += rankedPeaks.NotFoundTotal;
        }

        public RankProbability GetRank(int rankNum)
        {
            if (rankNum >= MaxRanks)
                rankNum = MaxRanks - 1;
            return (rankNum > 0 ? _rankTable[rankNum-1] : NotFound);
        }

        public Probability<IonType>[,] GetProbabilities()
        {
            var ranks = new Probability<IonType>[_rankTable.Count + 1, _ionTypes.Length];
            for (int i = 0; i < _rankTable.Count; i++)
            {
                for (int j = 0; j < _ionTypes.Length; j++)
                {
                    var ionType = _ionTypes[j];
                    ranks[i, j] = new Probability<IonType>(ionType, _rankTable[i].IonFrequencies[ionType], _rankTable[i].RankCount);
                }
            }
            for (int j = 0; j < _ionTypes.Length; j++)
            {
                var ionType = _ionTypes[j];
                ranks[ranks.GetLength(0)-1, j] = new Probability<IonType>(ionType, NotFound.IonFrequencies[ionType], NotFound.RankCount);
            }
            return ranks;
        }

        public IonType[] GetBinEdges()
        {
            return _ionTypes;
        }
    }
}
