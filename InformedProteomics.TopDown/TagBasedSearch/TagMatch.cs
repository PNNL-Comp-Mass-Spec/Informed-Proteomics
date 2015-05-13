namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class TagMatch
    {
        public TagMatch(int startIndex, int endIndex, int matchedTagLength, 
            int charge, double nTermScore, double cTermScore, double mass, string modifications)
        {
            StartIndex = startIndex;
            EndIndex = endIndex;
            MatchedTagLength = matchedTagLength;
            Charge = charge;
            NTermScore = nTermScore;
            CTermScore = cTermScore;
            Mass = mass;
            Modifications = modifications;
        }

        public int StartIndex { get; private set; } // Inclusive
        public int EndIndex { get; private set; }   // Exclusive
        public int MatchedTagLength { get; private set; }
        public int Charge { get; private set; }
        public double NTermScore { get; private set; }
        public double CTermScore { get; private set; }

        public double Score
        {
            get { return NTermScore + CTermScore + MatchedTagLength*2; }
        }

        public double Mass { get; private set; }
        public string Modifications { get; private set; }
    }
}
