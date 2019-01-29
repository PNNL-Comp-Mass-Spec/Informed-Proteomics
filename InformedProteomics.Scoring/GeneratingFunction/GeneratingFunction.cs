using System;
using System.Collections.Generic;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class GeneratingFunction
    {
        public GeneratingFunction(IScoringGraph graph)
        {
            _graph = graph;
            /*
            _maxAchievableScore = new int[_graph.GetNumNodes()];
            for (var i = _graph.GetNumNodes() - 1; i >= 0; i--)
            {
                var s = _graph.GetNodeScore(i);
                if (s > 0) _maxAchievableScore[i] = s;
                if (i < _graph.GetNumNodes() - 1) _maxAchievableScore[i] += _maxAchievableScore[i + 1];
            }*/
        }

        public double GetSpectralEValue(double score)
        {
            return _scoreDistribution.GetSpectralEValue(score);
        }

        public double GetSpectralEValue(int score)
        {
            return _scoreDistribution.GetSpectralEValue(score);
        }

        public ScoreDistribution GetScoreDistribution()
        {
            return _scoreDistribution;
        }
        //private int _targetScore;
        //private readonly int[] _maxAchievableScore;

        public void ComputeGeneratingFunction()//int targetScore = 0)
        {
            //_targetScore = targetScore;
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
            // ReSharper disable CommentTypo
            // TODO: adjusting the distribution depending on neighboring amino acid (e.g. R.PEPTIDEK. vs A.PEPTIDEK)
            // ReSharper restore CommentTypo
            _scoreDistribution = gfTable[gfTable.Length - 1];
        }

        public int MaximumScore => _scoreDistribution.MaxScore;

        private readonly IScoringGraph _graph;
        private ScoreDistribution _scoreDistribution;

        private ScoreDistribution GetScoreDistribution(int nodeIndex, IReadOnlyList<ScoreDistribution> gfTable)
        {
            var curNodeScore = _graph.GetNodeScore(nodeIndex);

            // determine MinScore and MaxScore
            var maxScore = (double) int.MinValue;
            var minScore = (double) int.MaxValue;

            var hasValidEdge = false;
            foreach (var edge in _graph.GetEdges(nodeIndex))
            {
                var prevNodeIndex = edge.PrevNodeIndex;
                var prevScoreDistribution = gfTable[prevNodeIndex];
                if (prevScoreDistribution == null) continue;
                var combinedScore = curNodeScore + _graph.GetEdgeScore(prevNodeIndex, nodeIndex);

                if (prevScoreDistribution.MaxScore + combinedScore > maxScore)
                {
                    maxScore = prevScoreDistribution.MaxScore + combinedScore;
                }
                if (prevScoreDistribution.MinScore + combinedScore < minScore)
                {
                    minScore = prevScoreDistribution.MinScore + combinedScore;
                }
                hasValidEdge = true;
            }

            if (!hasValidEdge) return null; //new ScoreDistribution(0, 1);

            // considering Target score for fast computing.....max 2 times faster but not useful at this point ScoreDistribution
            //minScore = Math.Max(_targetScore - _maxAchievableScore[nodeIndex], minScore);

            // Compute scoring distribution for the current node
            var scoringDistribution = new ScoreDistribution((int)Math.Floor(minScore), (int)Math.Ceiling(maxScore));

            foreach (var edge in _graph.GetEdges(nodeIndex))
            {
                var prevScoreDistribution = gfTable[edge.PrevNodeIndex];
                if (prevScoreDistribution == null) continue;
                var combinedScore = curNodeScore + _graph.GetEdgeScore(edge.PrevNodeIndex, nodeIndex);
                scoringDistribution.AddEValueDist(prevScoreDistribution, (int)Math.Round(combinedScore), edge.Weight);
            }
            return scoringDistribution;
        }
    }
    /*
    public class GeneratingFunction2
    {
        private readonly double[][] _eValueTable;
        private readonly int[] _minScoreAtNode;
        private readonly int[] _maxScoreAtNode;
        private bool _tableInit;
        private IScoringGraph _graph;

        public GeneratingFunction2(int maxPossibleNodes)
        {
            _eValueTable = new double[maxPossibleNodes][];
            _maxScoreAtNode = new int[maxPossibleNodes];
            _minScoreAtNode = new int[maxPossibleNodes];
            _eValueTable[0] = new double[1] {1}; // initialization for DP

            for (var i = 1; i < _eValueTable.Length; i++) _eValueTable[i] = new double[50];

            _tableInit = true;
        }

        public int MaximumScore
        {
            get
            {
                var lasNodeIndex = _graph.GetNumNodes() - 1;
                return _maxScoreAtNode[lasNodeIndex];
            }
        }

        public double[] GetSpectralEValueDistribution()
        {
            if (_graph == null) return null;

            var lasNodeIndex = _graph.GetNumNodes() - 1;
            var lastRow = _eValueTable[lasNodeIndex];
            var maxScore = _maxScoreAtNode[lasNodeIndex];
            var minScore = _minScoreAtNode[lasNodeIndex];

            var eValueDist = new double[maxScore - minScore + 1];
            Array.Copy(lastRow, eValueDist, eValueDist.Length);

            return eValueDist;
        }

        public double GetSpectralEValue(double[] eValueDist, int score)
        {
            if (_graph == null) return double.MaxValue;

            var spectralEValue = 0d;
            for (var col = score; col < eValueDist.Length; col++)
            {
                spectralEValue += eValueDist[col];
            }
            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        public double GetSpectralEValue(int score)
        {
            if (_graph == null) return double.MaxValue;

            var spectralEValue = 0d;
            var lasNodeIndex = _graph.GetNumNodes() - 1;
            var lastRow = _eValueTable[lasNodeIndex];
            for (var col = score; col <= _maxScoreAtNode[lasNodeIndex]; col++)
            {
                spectralEValue += lastRow[col];
            }

            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        public void ComputeGeneratingFunction(IScoringGraph graph)
        {
            if (graph == null) return;

            _graph = graph;

            InitTable();
            for (var curNodeIndex = 1; curNodeIndex < _graph.GetNumNodes(); curNodeIndex++)
            {
                var curNodeScore = _graph.GetNodeScore(curNodeIndex);

                // find min/max scores of current node
                var minScore = int.MaxValue;
                var maxScore = int.MinValue;
                var hasValidEdge = false;
                foreach (var edge in _graph.GetEdges(curNodeIndex))
                {
                    var prevNodeIndex = edge.PrevNodeIndex;
                    var combinedScore = curNodeScore + edge.Score;
                    maxScore = Math.Max(_maxScoreAtNode[prevNodeIndex] + combinedScore, maxScore);
                    minScore = Math.Min(_minScoreAtNode[prevNodeIndex] + combinedScore, minScore);
                    hasValidEdge = true;
                }

                if (!hasValidEdge) continue;

                _minScoreAtNode[curNodeIndex] = minScore;
                _maxScoreAtNode[curNodeIndex] = maxScore;

                if (_eValueTable[curNodeIndex].Length < maxScore - minScore + 1) _eValueTable[curNodeIndex] = new double[maxScore - minScore + 1];

                var curNodeEValues = _eValueTable[curNodeIndex];
                foreach (var edge in _graph.GetEdges(curNodeIndex))
                {
                    var prevNodeIndex = edge.PrevNodeIndex;
                    var combinedScore = curNodeScore + edge.Score;

                    var prevNodeEValues = _eValueTable[prevNodeIndex];
                    var prevNodeMinScore = _minScoreAtNode[prevNodeIndex];

                    var effectiveMinScoreAtNode = Math.Max(prevNodeMinScore + combinedScore, _minScoreAtNode[curNodeIndex]);
                    var effectiveMaxScoreAtNode = _maxScoreAtNode[prevNodeIndex] + combinedScore;

                    for (var score = effectiveMinScoreAtNode; score <= effectiveMaxScoreAtNode; score++)
                    {
                        curNodeEValues[score - minScore] += prevNodeEValues[score - combinedScore - prevNodeMinScore] * edge.Weight;
                    }
                }
            }

            _tableInit = false;
        }

        private void InitTable()
        {
            if (_tableInit) return;
            //Array.Clear(_eValueTable, 1, _eValueTable.Length - 1);
            for (var i = 1; i < _eValueTable.Length; i++) Array.Clear(_eValueTable[i], 0, _eValueTable[i].Length);

            Array.Clear(_minScoreAtNode, 0, _minScoreAtNode.Length);
            Array.Clear(_maxScoreAtNode, 0, _maxScoreAtNode.Length);
            _tableInit = true;
        }
    }*/
}