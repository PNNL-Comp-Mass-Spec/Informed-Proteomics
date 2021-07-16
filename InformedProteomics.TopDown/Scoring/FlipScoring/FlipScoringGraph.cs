using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Class for scoring an amino acid graph using the FLIP scoring model.
    /// </summary>
    public class FlipScoringGraph : IScoringGraph
    {
        /// <summary>
        /// The mass bins for constructing the nodes of the graph.
        /// </summary>
        private readonly IMassBinning massBins;

        /// <summary>
        /// Node scores for nTerminal ions.
        /// </summary>
        private readonly double?[] nTerminalNodes;

        /// <summary>
        /// Node scores for cTerminal ions.
        /// </summary>
        private readonly double?[] cTerminalNodes;

        /// <summary>
        /// Mapping between sink node -> graph edge.
        /// </summary>
        private readonly Dictionary<int, List<FlipScoringGraphEdge>> edges;

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScoringGraph" />.
        /// </summary>
        /// <param name="massBins">The mass bins for constructing the nodes of the graph.</param>
        /// <param name="nTerminalNodes">Node scores for nTerminal ions.</param>
        /// <param name="cTerminalNodes">Node scores for cTerminal ions.</param>
        /// <param name="edges">Mapping between sink node -> graph edge.</param>
        public FlipScoringGraph(
                                IMassBinning massBins,
                                double?[] nTerminalNodes,
                                double?[] cTerminalNodes,
                                Dictionary<int, List<FlipScoringGraphEdge>> edges)
        {
            this.massBins = massBins;
            this.nTerminalNodes = nTerminalNodes;
            this.cTerminalNodes = cTerminalNodes;
            this.edges = edges;
        }

        /// <summary>
        /// Gets the node score for a particular sink node index.
        /// Node score is the ion score WITHOUT the mass error score or amino acid probability score.
        /// </summary>
        /// <param name="nodeIndex">The index (nominal mass) of the sink node.</param>
        /// <returns>The node score of the sink node.</returns>
        public double GetNodeScore(int nodeIndex)
        {
            return nTerminalNodes[nodeIndex] ?? 0 + cTerminalNodes[nodeIndex] ?? 0;
        }

        /// <summary>
        /// Get all of the in edges for a particular sink node.
        /// </summary>
        /// <param name="nodeIndex">The index (nominal mass) of the sink node.</param>
        /// <returns><see cref="List{T}" /> of <see ref="FlipScoringGraphEdge" />.</returns>
        public IEnumerable<IScoringGraphEdge> GetEdges(int nodeIndex)
        {
            return edges[nodeIndex];
        }

        /// <summary>
        /// Get edge score
        /// </summary>
        /// <param name="sourceMassBin"></param>
        /// <param name="sinkMassBin"></param>
        /// <returns>Score for the edge between the two bins</returns>
        public double GetEdgeScore(int sourceMassBin, int sinkMassBin)
        {
            if (GetEdges(sinkMassBin).First(e => e.PrevNodeIndex == sourceMassBin) is not FlipScoringGraphEdge edge)
            {
                return 0.0;
            }

            return edge.GetScore();
        }

        /// <summary>
        /// Gets the score for an entire sequence against the scoring graph
        /// </summary>
        /// <param name="sequence">The sequence to score.</param>
        /// <returns>Score</returns>
        public double ScoreSequence(Sequence sequence)
        {
            var cleavages = sequence.GetInternalCleavages();
            var score = 0.0;
            int prevNTermBin = 0, prevCTermBin = 0;
            foreach (var cleavage in cleavages)
            {
                var nTermBin = massBins.GetBinNumber(cleavage.PrefixComposition.Mass);
                var cTermBin = massBins.GetBinNumber(cleavage.SuffixComposition.Mass);

                score += GetNodeScore(nTermBin) + GetNodeScore(cTermBin);
                score += GetEdgeScore(prevNTermBin, nTermBin) + GetEdgeScore(prevCTermBin, cTermBin);

                prevNTermBin = nTermBin;
                prevCTermBin = cTermBin;
            }

            return score;
        }

        /// <summary>
        /// Gets the total number of nodes (mass bins) in the scoring graph.
        /// </summary>
        /// <returns>The total number of nodes (mass bins) in the scoring graph.</returns>
        public int GetNumNodes()
        {
            return nTerminalNodes.Length;
        }
    }
}
