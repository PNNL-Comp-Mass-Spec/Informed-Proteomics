using System.Collections.Generic;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public interface IScoringGraph
    {
        double GetNodeScore(int nodeIndex);
        IEnumerable<IScoringGraphEdge> GetEdges(int nodeIndex);

        double GetEdgeScore(int nodeIndex1, int nodeIndex2);
        int GetNumNodes();
    }
}
