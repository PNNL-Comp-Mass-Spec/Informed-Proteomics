using System.Collections.Generic;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public interface IScoringGraph
    {
        int GetNodeScore(int nodeIndex);
        IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex);
        int GetNumNodes();
    }
}
