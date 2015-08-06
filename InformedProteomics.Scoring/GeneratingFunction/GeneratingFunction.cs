using System;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class GeneratingFunction
    {
        public GeneratingFunction(int maxPossibleNodes)
        {
            _eValueTable = new double[maxPossibleNodes][];
            _maxScoreAtNode = new int[maxPossibleNodes];
            _minScoreAtNode = new int[maxPossibleNodes];
            _eValueTable[0] = new double[1] { 1 }; // initialization for DP
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
            var lasNodeIndex = _graph.GetNumNodes() - 1;
            var lastRow = _eValueTable[lasNodeIndex];
            var maxScore = _maxScoreAtNode[lasNodeIndex];
            var evalueDist = new double[maxScore + 1];
            Array.Copy(lastRow, evalueDist, evalueDist.Length);

            return evalueDist;
        }

        public double GetSpectralEValue(double[] evalueDist, int score)
        {
            var spectralEValue = 0d;
            for (var col = score; col < evalueDist.Length; col++)
            {
                spectralEValue += evalueDist[col];
            }
            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        public double GetSpectralEValue(int score)
        {
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

                _minScoreAtNode[curNodeIndex] = (hasValidEdge) ? minScore : 0;
                _maxScoreAtNode[curNodeIndex] = (hasValidEdge) ? maxScore : 1;

                _eValueTable[curNodeIndex] = new double[(hasValidEdge) ? maxScore + 1 : 2];
                var curNodeEvalues = _eValueTable[curNodeIndex];


                foreach (var edge in _graph.GetEdges(curNodeIndex))
                {
                    var prevNodeIndex = edge.PrevNodeIndex;
                    var combinedScore = curNodeScore + edge.Score;

                    var prevNodeEvalues = _eValueTable[prevNodeIndex];
                    var effectiveMinScoreAtNode = Math.Max(_minScoreAtNode[prevNodeIndex] + combinedScore, _minScoreAtNode[curNodeIndex]);
                    var effectiveMaxScoreAtNode = _maxScoreAtNode[prevNodeIndex] + combinedScore;

                    for (var score = effectiveMinScoreAtNode; score <= effectiveMaxScoreAtNode; score++)
                    {
                        curNodeEvalues[score] += prevNodeEvalues[score - combinedScore] * edge.Weight;
                    }
                }
            }
            
            _tableInit = false;
        }
      
        private void InitTable()
        {
            if (_tableInit) return;
            Array.Clear(_eValueTable, 1, _eValueTable.Length - 1);
            Array.Clear(_minScoreAtNode, 0, _minScoreAtNode.Length);
            Array.Clear(_maxScoreAtNode, 0, _maxScoreAtNode.Length);
            _tableInit = true;
        }
        
        private readonly double[][] _eValueTable;
        private readonly int[] _minScoreAtNode;
        private readonly int[] _maxScoreAtNode;
        private bool _tableInit;
        private IScoringGraph _graph;
    }
}