namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoringGraphEdge
    {
        public ScoringGraphEdge(int prevNodeIndex, int score, double weight)
        {
            PrevNodeIndex = prevNodeIndex;
            //Score = score;
            Weight = weight;
        }

        public int PrevNodeIndex { get; private set; }
        //public int Score { get; private set; }
        public int Score { get { return 0; } }
        public double Weight { get; private set; }
    }
}
