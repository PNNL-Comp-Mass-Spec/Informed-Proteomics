namespace InformedProteomics.Scoring.GeneratingFunction
{
    public interface IScoringGraphEdge
    {
        int PrevNodeIndex { get; }
        double Weight { get; }
    }
}
