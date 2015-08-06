using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    /*
    public class ScoringGraph : IScoringGraph
    {
        
        public ScoringGraph(double[] massList, int[] peakScore, double proteinMass, AminoAcid[] aminoArray, double[] aminoAcidProb, ProteinMassComparerWithBinning comparer)
        {
            var maxBinIndex = comparer.GetBinNumber(proteinMass);
            _aminoArray = aminoArray;


            var stopwatch = Stopwatch.StartNew();

            //var graph = new ScoringGraph(massList, peakScores, proteinMass, aminoArray, aminoProb, comparer);
            
            // node generation
            _scoreNodes = new ScoringGraphNode[maxBinIndex + 1];
            for (var j = 0; j < massList.Length; j++)
            {
                var i = comparer.GetBinNumber(massList[j]);
                if (i < 0 || i > maxBinIndex) continue; // ignore mass peak
                _scoreNodes[i] = new ScoringGraphNode(peakScore[j]);
            }

            for (var i = 0; i <= maxBinIndex; i++)
            {
                if (_scoreNodes[i] == null) _scoreNodes[i] = new ScoringGraphNode();
            }
            Console.WriteLine(@"node generation elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

            stopwatch.Restart();
            // edge generation
            for (var i = 0; i <= maxBinIndex; i++)
            {
                var nodeMass = comparer.GetMass(i);
                for (var a = 0; a < aminoArray.Length; a++)
                {
                    var j = comparer.GetBinNumber(nodeMass + aminoArray[a].Mass);
                    if (j < 0 || j > maxBinIndex) continue;
                    var edge = new ScoringGraphEdge(i, EdgeScore, aminoAcidProb[a]);
                    _scoreNodes[j].AddEdge(edge);
                }
            }
            Console.WriteLine(@"node generation elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
        }

        public ScoringGraph(double[] massList, int[] peakScore, double proteinMass, AminoAcid[] aminoArray, double[] aminoAcidProb)
        {
            var maxBinIndex = Constants.GetBinNumHighPrecision(proteinMass);
            _aminoArray = aminoArray;

            // node generation
            _scoreNodes = new ScoringGraphNode[maxBinIndex + 1];
            for (var j = 0; j < massList.Length; j++)
            {
                var i = Constants.GetBinNumHighPrecision(massList[j]);
                if (i < 0 || i > maxBinIndex) continue; // ignore mass peak
                _scoreNodes[i] = new ScoringGraphNode(peakScore[j]);
            }

            for (var i = 0; i <= maxBinIndex; i++)
            {
                if (_scoreNodes[i] == null) _scoreNodes[i] = new ScoringGraphNode();
            }

            // edge generation
            for (var i = 0; i <= maxBinIndex; i++)
            {
                var nodeMass = i/Constants.RescalingConstantHighPrecision;
                for (var a = 0; a < aminoArray.Length; a++)
                {
                    var j = Constants.GetBinNumHighPrecision(nodeMass + aminoArray[a].Mass);
                    if (j < 0 || j > maxBinIndex) continue;
                    var edge = new ScoringGraphEdge(i, EdgeScore, aminoAcidProb[a]);
                    _scoreNodes[j].AddEdge(edge);
                }
            }
        }

        public int GetNodeScore(int nodeIndex)
        {
            return _scoreNodes[nodeIndex].Score;
        }

        public IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex)
        {
            return _scoreNodes[nodeIndex].Edges;
        }

        public int GetNumNodes()
        {
            return _scoreNodes.Length;
        }

        private readonly AminoAcid[] _aminoArray;
        private const int EdgeScore = 0;
        private readonly ScoringGraphNode[] _scoreNodes;
    }
    */
}
