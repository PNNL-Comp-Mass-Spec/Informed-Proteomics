using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.GeneratingFunction
{

    public class GeneratingFunction
    {
        public GeneratingFunction(IScoringGraph graph)
        {
            _graph = graph;
        }

        public void ComputeGeneratingFunction()
        {
            throw new NotImplementedException();
        }

        public double GetSpectralProbability(int score)
        {
            throw new NotImplementedException();
        }
        
        public double GetSpectralEValue(int score)
        {
            throw new NotImplementedException();
        }

        private readonly IScoringGraph _graph;
    }





/*
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

        private readonly IScoringGraph _graph;
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
                if (prevScoreDistribution.MinScore + combinedScore > minScore)
                {
                    minScore = prevScoreDistribution.MinScore + combinedScore;
                }
                validEdges.Add(edge);
            }

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
*/
}
