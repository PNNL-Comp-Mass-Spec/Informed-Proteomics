using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class TagMatch
    {
        public TagMatch(int startIndex, int endIndex, int matchedTagLength,
            int charge, double nTermScore, double cTermScore, double mass, ModificationCombination modifications, string modificationText)
        {
            StartIndex = startIndex;
            EndIndex = endIndex;
            MatchedTagLength = matchedTagLength;
            Charge = charge;
            NTermScore = nTermScore;
            CTermScore = cTermScore;
            Mass = mass;
            Modifications = modifications;
            ModificationText = modificationText;
        }

        public int StartIndex { get; } // Inclusive
        public int EndIndex { get; }   // Exclusive
        public int MatchedTagLength { get; }
        public int Charge { get; }
        public double NTermScore { get; }
        public double CTermScore { get; }

        public double Score { get; internal set; }
        /*
        public double Score
        {
            get { return NTermScore + CTermScore + MatchedTagLength*2; }
        }*/

        public double Mass { get; }
        public ModificationCombination Modifications { get; }
        public string ModificationText { get; }
    }
}
