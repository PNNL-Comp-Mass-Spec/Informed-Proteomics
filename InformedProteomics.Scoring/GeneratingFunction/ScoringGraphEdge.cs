namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoringGraphEdge : IScoringGraphEdge
    {
        public ScoringGraphEdge(int prevNodeIndex)
        {
            PrevNodeIndex = prevNodeIndex;
            //Score = score;
        }

        public int PrevNodeIndex { get; private set; }
        //public int Score { get; private set; }
        public double Weight { get { return 0.05; } }
    }
}
