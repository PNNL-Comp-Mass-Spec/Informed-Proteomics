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
            if (_graph == null) return null;
            
            var lasNodeIndex = _graph.GetNumNodes() - 1;
            var lastRow = _eValueTable[lasNodeIndex];
            var maxScore = _maxScoreAtNode[lasNodeIndex];
            var evalueDist = new double[maxScore + 1];
            Array.Copy(lastRow, evalueDist, evalueDist.Length);

            return evalueDist;
        }

        public double GetSpectralEValue(double[] evalueDist, int score)
        {
            if (_graph == null) return double.MaxValue;

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
            _graph = graph;
            
            if (graph == null) return;

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

    /* 
     * Sangtae's implementation
    public class GeneratingFunction
    {
        public GeneratingFunction(IScoringGraph graph)
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
    */
}