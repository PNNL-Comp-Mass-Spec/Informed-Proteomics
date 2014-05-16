namespace InformedProteomics.BottomUp.Scoring
{
    public class IcBottomUpScores
    {
        public IcBottomUpScores(double score, string modifications)
        {
            Score = score;
            Modifications = modifications;
        }

        public double Score { get; private set; }
        public string Modifications { get; private set; }
    }
}
