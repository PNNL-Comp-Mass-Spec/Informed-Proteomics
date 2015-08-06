using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Scoring.GeneratingFunction;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class GeneratingFunctionOld
    {
        public GeneratingFunctionOld(IScoringGraph graph)
        {
            _graph = graph;
        }        
        
        public double GetSpectralEValue(int score)
        {
            return _scoreDistribution.GetSpectralEValue(score);
        }

        public void ComputeGeneratingFunction()
        {
            var gfTable = new ScoreDistribution[_graph.GetNumNodes()];
            // Source
            var sourceDist = new ScoreDistribution(0, 1);
            sourceDist.SetEValue(0, 1);
            gfTable[0] = sourceDist;

            // All the other nodes
            for (var nodeIndex = 1; nodeIndex < _graph.GetNumNodes(); nodeIndex++)
            {
                gfTable[nodeIndex] = GetScoreDistribution(nodeIndex, gfTable);
            }

            // Sink
            // TODO: adjusting the distribution depending on neighboring amino acid (e.g. R.PEPTIDEK. vs A.PEPTIDEK)

            _scoreDistribution = gfTable[gfTable.Length - 1];
        }

        
        private IScoringGraph _graph;
        private ScoreDistribution _scoreDistribution;

        private ScoreDistribution GetScoreDistribution(int nodeIndex, ScoreDistribution[] gfTable)
        {
            var curNodeScore = _graph.GetNodeScore(nodeIndex);

            // determine MinScore and MaxScore
            var maxScore = int.MinValue;
            var minScore = int.MaxValue;

            var validEdges = new List<ScoringGraphEdge>();

            foreach (var edge in _graph.GetEdges(nodeIndex))
            {
                var prevNodeIndex = edge.PrevNodeIndex;
                var prevScoreDistribution = gfTable[prevNodeIndex];
                if (prevScoreDistribution == null) continue;
                var combinedScore = curNodeScore + edge.Score;
                if (prevScoreDistribution.MaxScore + combinedScore > maxScore)
                {
                    maxScore = prevScoreDistribution.MaxScore + combinedScore;
                }
                if (prevScoreDistribution.MinScore + combinedScore < minScore)
                {
                    minScore = prevScoreDistribution.MinScore + combinedScore;
                }
                validEdges.Add(edge);
            }

            if (validEdges.Count < 1) return new ScoreDistribution(0, 1);

            //if (nodeIndex % 1000 == 0) Console.WriteLine("{0}\t{1}\t{2}",nodeIndex, minScore, maxScore);

            // Compute scoring distribution for the current node
            var scoringDistribution = new ScoreDistribution(minScore, maxScore);
            foreach (var edge in validEdges)
            {
                var prevScoreDistribution = gfTable[edge.PrevNodeIndex];
                if (prevScoreDistribution == null) continue;
                var combinedScore = curNodeScore + edge.Score;
                scoringDistribution.AddEValueDist(prevScoreDistribution, combinedScore, edge.Weight);
            }
            return scoringDistribution;
        }
    }
}
