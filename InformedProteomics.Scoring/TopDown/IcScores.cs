namespace InformedProteomics.Scoring.TopDown
{
    public class IcScores
    {
        public IcScores(int nMatchedFragments, double score, string modifications)
        {
            NumMatchedFrags = nMatchedFragments;
            Score = score;
            Modifications = modifications;
        }

        public int NumMatchedFrags { get; }
        public double Score { get; } // this score is used to calculate p-value by generating function

        public string Modifications { get; }

        public override string ToString()
        {
            return string.Join("\t",
                new[]
                {
                    NumMatchedFrags, Score,
                });
        }

        public static string GetScoreNames()
        {
            return "#MatchedFragments\tScore";
        }
    }
}
