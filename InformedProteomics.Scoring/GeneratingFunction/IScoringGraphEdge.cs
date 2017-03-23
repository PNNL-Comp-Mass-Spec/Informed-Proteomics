using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public interface IScoringGraphEdge
    {
        int PrevNodeIndex { get; }
        double Weight { get; }
    }
}
