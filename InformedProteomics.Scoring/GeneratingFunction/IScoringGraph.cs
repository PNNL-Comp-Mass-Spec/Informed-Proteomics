using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public interface IScoringGraph
    {
        int GetNodeScore(int nodeIndex);
        IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex);
        int GetNumNodes();
    }
}
