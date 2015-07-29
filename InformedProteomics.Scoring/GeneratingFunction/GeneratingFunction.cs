using System;
using System.Collections.Generic;
using System.Globalization;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class GeneratingFunction
    {
        public GeneratingFunction(int maxPossibleNodes, int maxPossibleScore)
        {
            _eValueTable = new double[maxPossibleScore][];
            for (var i = 0; i < maxPossibleScore; i++)
            {
                _eValueTable[i] = new double[maxPossibleNodes];
            }

            _eValueTable[0][0] = 1;
            _tableInit = true;
        }

        public void ComputeGeneratingFunction(IScoringGraph graph)
        {
            _graph = graph;
            InitTable();

            for (var i = 0; i < _eValueTable.Length; i++)
            {
                for (var j = 1; j < _graph.GetNumNodes(); j++)
                {
                    var evalue = 0;

                    _eValueTable[i][j] = evalue;
                }

                if (_eValueTable[i][_graph.GetNumNodes() - 1] < double.Epsilon)
                {
                    _maxScore = i;
                    break;
                }
            }

            _tableInit = false;
        }

        public double GetSpectralEValue(int score)
        {
            var spectralEValue = 0.0;
            for (var i =0; i <= _maxScore; i++)
            {
                spectralEValue += _eValueTable[i][_graph.GetNumNodes() - 1];
            }
            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        private void InitTable()
        {
            if (_tableInit) return;
            Array.Clear(_eValueTable, 0, _eValueTable.Length);
            for (var i = 0; i < _eValueTable.Length; i++) Array.Clear(_eValueTable[i], 0, _eValueTable[i].Length);
            _eValueTable[0][0] = 1;
            _tableInit = true;
        }


        private readonly double[][] _eValueTable;
        private bool _tableInit;
        
        private IScoringGraph _graph;
        private int _maxScore;

        /*
        public GeneratingFunction(IScoringGraph graph)
        {
            _graph = graph;
            _eValueTable = new List<double[]>();
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
}
