using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.TopDown.Execution
{
    public class ProteoformSpectrumMatchContainer
    {
        public ProteoformSpectrumMatchContainer(FastaDatabase database, int[] ms2ScanVector, int maxModifications, int maxNumMatchesPerSpectrum, int minScore = 4)
        {
            Database = database;
            NumMatchesPerSpectrum = maxNumMatchesPerSpectrum;
            _scoreCutoff = minScore;
            Ms2ScanVector = ms2ScanVector;
            _ms2ScanToIndexMap = new int[ms2ScanVector.Last() + 1];
            for(var i = 0; i < ms2ScanVector.Length; i++)
            {
                var scanNum = ms2ScanVector[i];
                _ms2ScanToIndexMap[scanNum] = i;
            }
            _matchedSet = new SortedSet<DatabaseSequenceSpectrumMatch>[maxModifications + 1][];
            for(var i = 0; i <= maxModifications; i++) _matchedSet[i] = new SortedSet<DatabaseSequenceSpectrumMatch>[ms2ScanVector.Length];

            _checkedOutScanNumbers = new List<int>();
        }

        public void AddMatch(DatabaseSequenceSpectrumMatch newMatch)
        {
            if (newMatch.Score < _scoreCutoff) return;
            var scanIndex = _ms2ScanToIndexMap[newMatch.ScanNum];
            var modIndex = (newMatch.Modifications == null) ? 0 : newMatch.Modifications.GetNumModifications();
            if (modIndex >= _matchedSet.Length) return;

            // thread safe
            lock (_matchedSet[modIndex])
            {
                if (_matchedSet[modIndex][scanIndex] == null)
                {
                    _matchedSet[modIndex][scanIndex] = new SortedSet<DatabaseSequenceSpectrumMatch> { newMatch };
                }
                else // already exists
                {
                    var existingMatches = _matchedSet[modIndex][scanIndex];
                    var maxScore = existingMatches.Max.Score;
                    if (existingMatches.Count < NumMatchesPerSpectrum && maxScore * ScoreRatioCutoff < newMatch.Score)
                    {
                        existingMatches.Add(newMatch);
                        existingMatches.RemoveWhere(mt => mt.Score < maxScore * ScoreRatioCutoff);
                    }
                    else
                    {
                        var minScore = existingMatches.Min.Score;
                        if (newMatch.Score > minScore)
                        {
                            existingMatches.Add(newMatch);
                            existingMatches.RemoveWhere(mt => mt.Score < maxScore * ScoreRatioCutoff);
                        }
                    }
                }
            }
        }

        public int GetNumberOfMatches()
        {
            return (from matchesWithSameNMods in _matchedSet where matchesWithSameNMods != null from matchesPerSpec in matchesWithSameNMods where matchesPerSpec != null select matchesPerSpec.Count).Sum();
        }

        public SortedSet<DatabaseSequenceSpectrumMatch>[] GetMatches(int nModifications)
        {
            return _matchedSet[nModifications];
        }

        private const double ScoreRatioCutoff = 0.7;
        public readonly int NumMatchesPerSpectrum;
        public readonly int[] Ms2ScanVector;
        private readonly int[] _ms2ScanToIndexMap;
        public readonly FastaDatabase Database;
        private readonly int _scoreCutoff;
        private readonly SortedSet<DatabaseSequenceSpectrumMatch>[][] _matchedSet;

        private readonly List<int> _checkedOutScanNumbers;

        //private readonly HashSet<int> checkOutFlagMs2Scans;
    }
}
