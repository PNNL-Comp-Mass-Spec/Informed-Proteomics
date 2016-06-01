using System.Collections.Generic;

namespace InformedProteomics.TopDown.Scoring
{
    public class SubScoreFactory
    {
        private readonly Dictionary<int, int> _missingXicScore; // charge, score
        private readonly Dictionary<int, Dictionary<int, int>> _xicCorrScore;    // charge, correlation raw score, score

        private void Read(string filePath)
        {

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
