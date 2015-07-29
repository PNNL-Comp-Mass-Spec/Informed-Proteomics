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
            _edges = null;
        }

        public void AddEdge(ScoringGraphEdge edge)
        {
            if (_edges == null) _edges = new List<ScoringGraphEdge>();
            _edges.Add(edge);
        }

        public IEnumerable<ScoringGraphEdge> Edges
        {
            get
            {
                for (var i = 0; _edges != null && i < _edges.Count; i++)
                {
                    yield return _edges[i];
                }
            }
        }

        public readonly int Score;
        private List<ScoringGraphEdge> _edges;
    }
}
