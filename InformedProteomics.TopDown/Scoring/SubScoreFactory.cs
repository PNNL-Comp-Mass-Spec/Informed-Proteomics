using System;
using System.Collections.Generic;

namespace InformedProteomics.TopDown.Scoring
{
    [Obsolete("This functionality of this class is not implemented")]
    public class SubScoreFactory
    {
        // charge, score
        private readonly Dictionary<int, int> _missingXicScore;

        // charge, correlation raw score, score
        private readonly Dictionary<int, Dictionary<int, int>> _xicCorrScore;

        /// <summary>
        /// Constructor
        /// </summary>
        public SubScoreFactory()
        {
            _missingXicScore = new Dictionary<int, int>();
            _xicCorrScore = new Dictionary<int, Dictionary<int, int>>();
        }

        public int GetMissingXicScore(int charge)
        {
            var score = _missingXicScore[charge];
            return score;
        }

        public int GetXicCorrScore(int charge, int correlationIntegerRawScore)
        {
            var t = _xicCorrScore[charge];
            if (t == null) return 0;
            return t[correlationIntegerRawScore];
        }
    }
}
