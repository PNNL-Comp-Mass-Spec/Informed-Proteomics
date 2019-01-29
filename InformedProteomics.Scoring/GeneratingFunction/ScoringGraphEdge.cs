namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoringGraphEdge : IScoringGraphEdge
    {
        public ScoringGraphEdge(int prevNodeIndex)
        {
            PrevNodeIndex = prevNodeIndex;
            //Score = score;
        }

        public int PrevNodeIndex { get; }
        //public int Score { get; private set; }
        public double Weight => 0.05;
    }
}
