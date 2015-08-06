using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoringGraphNode
    {
        public ScoringGraphNode(int nodeScore = 0)
        {
            Score = nodeScore;
            _edgeList = null;
        }
        
        public void AddEdge(ScoringGraphEdge edge)
        {
            if (_edgeList == null) _edgeList = new List<ScoringGraphEdge>();
            _edgeList.Add(edge);
        }

        public IEnumerable<ScoringGraphEdge> Edges
        {
            get
            {
                for (var i = 0; _edgeList != null && i < _edgeList.Count; i++)
                {
                    yield return _edgeList[i];
                }
            }
        }
        
        public readonly int Score;
        private List<ScoringGraphEdge> _edgeList;
    }
}
