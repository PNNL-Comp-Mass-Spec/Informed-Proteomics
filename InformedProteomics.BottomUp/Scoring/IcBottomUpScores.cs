namespace InformedProteomics.BottomUp.Scoring
{
    public class IcBottomUpScores
    {
        public IcBottomUpScores(double score, string modifications)
        {
            Score = score;
            Modifications = modifications;
        }

        public double Score { get; }
        public string Modifications { get; }

        /// <inheritdoc />
        public override string ToString()
        {
            return Score+"\t"+Modifications;
        }
    }
}
