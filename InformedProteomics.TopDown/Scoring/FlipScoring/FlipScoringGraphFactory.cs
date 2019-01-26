using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Creates a scoring graph using the FLIP scoring model.
    /// </summary>
    public class FlipScoringGraphFactory
    {
        /// <summary>
        /// The mass bins for constructing the nodes of the graph.
        /// </summary>
        private readonly IMassBinning massBins;

        /// <summary>
        /// The pre-calculated edges for the scoring graph.
        /// </summary>
        private readonly List<FlipScoringGraphEdge> edges;

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScoringGraphFactory" />.
        /// </summary>
        /// <param name="massBins">The mass bins for constructing the nodes of the graph.</param>
        /// <param name="aminoAcidSet">Amino acid set to build the graph edges from.</param>
        /// <param name="aminoAcidProbabilities">The amino acid probabilities.</param>
        public FlipScoringGraphFactory(IMassBinning massBins, AminoAcidSet aminoAcidSet, Dictionary<char, double> aminoAcidProbabilities)
        {
            this.massBins = massBins;
            this.edges = this.InitEdges(aminoAcidSet, aminoAcidProbabilities);
        }

        /// <summary>
        /// Create a scoring graph that uses the FLIP scoring model.
        /// </summary>
        /// <param name="flipScorer">The FLIP scoring model to use.</param>
        /// <param name="proteinMass">The maximum mass of the scoring graph.</param>
        /// <returns>The scoring graph that uses the FLIP scoring model.</returns>
        public FlipScoringGraph GetScoringGraph(FlipScorer<DeconvolutedSpectrum> flipScorer, double proteinMass)
        {
            var proteinMassBin = this.massBins.GetBinNumber(proteinMass);

            // Precalculate scores for all nodes
            var nodeScores = this.InitNodeScores(flipScorer, proteinMass);

            // Filter the edges for the scoring graph.
            var graphEdges = this.edges.Where(edge => edge.SinkNodeIndex <= proteinMassBin)
                                       .Select(edge => new FlipScoringGraphEdge(edge.PrevNodeIndex, edge.SinkNodeIndex, edge.Weight, edge.Label, flipScorer))
                                       .GroupBy(edge => edge.SinkNodeIndex)
                                       .ToDictionary(edgeGroup => edgeGroup.Key, edgeGroup => edgeGroup.ToList());

            return new FlipScoringGraph(this.massBins, nodeScores[0], nodeScores[1], graphEdges);
        }

        /// <summary>
        /// Precompute edges for the scoring graph.
        /// </summary>
        /// <param name="aminoAcidSet">Amino acid set to build the graph edges from.</param>
        /// <param name="aminoAcidProbabilities">The amino acid probabilities.</param>
        /// <returns>A list of all scoring graph edges.</returns>
        private List<FlipScoringGraphEdge> InitEdges(AminoAcidSet aminoAcidSet, Dictionary<char, double> aminoAcidProbabilities)
        {
            var adjList = new LinkedList<FlipScoringGraphEdge>[this.massBins.NumberOfBins];
            for (var i = 0; i < this.massBins.NumberOfBins; i++) adjList[i] = new LinkedList<FlipScoringGraphEdge>();

            var terminalModifications = FilteredProteinMassBinning.GetTerminalModifications(aminoAcidSet);
            var aminoAcidArray = FilteredProteinMassBinning.GetExtendedAminoAcidArray(aminoAcidSet);

            for (var i = 0; i < this.massBins.NumberOfBins; i++)
            {
                var mi = this.massBins.GetMass(i);
                var fineNodeMass = mi;

                foreach (var aa in aminoAcidArray)
                {
                    var j = this.massBins.GetBinNumber(fineNodeMass + aa.Mass);
                    if (j < 0 || j >= this.massBins.NumberOfBins) continue;
                    var aaWeight = aminoAcidProbabilities.ContainsKey(aa.Residue) ? Math.Log10(aminoAcidProbabilities[aa.Residue]) : 0;
                    adjList[j].AddLast(new FlipScoringGraphEdge(i, j, aaWeight, aa, null));

                    if (i == 0 && !(aa is ModifiedAminoAcid))
                    {
                        foreach (var terminalMod in terminalModifications)
                        {
                            var modifiedAa = new ModifiedAminoAcid(aa, terminalMod);
                            j = this.massBins.GetBinNumber(fineNodeMass + modifiedAa.Mass);
                            if (j < 0 || j >= this.massBins.NumberOfBins) continue;
                            adjList[j].AddLast(new FlipScoringGraphEdge(i, j, aaWeight, modifiedAa, null));
                        }
                    }
                }
            }

            return adjList.SelectMany(edge => edge).ToList();
        }

        /// <summary>
        /// Get the scores of the nodes in the GeneratingFunction scoring graph.
        /// </summary>
        /// <param name="scorer">The FLIP scoring model to use.</param>
        /// <param name="proteinMass">The precursor mass.</param>
        private double?[][] InitNodeScores(FlipScorer<DeconvolutedSpectrum> scorer, double proteinMass)
        {
            var numNodes = this.massBins.GetBinNumber(proteinMass) + 1;
            var nodeScores = new double?[2][];
            var nTerminalFragScores = nodeScores[0] = new double?[numNodes];
            var cTerminalFragScores = nodeScores[1] = new double?[numNodes];

            var deconvPeaks = (DeconvolutedPeak[]) scorer.ProductSpectrum.Peaks;

            // transform all peaks into being nterminal ions
            foreach (var peak in deconvPeaks)
            {
                foreach (var ionType in scorer.SelectedIonTypes)
                {
                    var scores = ionType.IsPrefix ? nTerminalFragScores : cTerminalFragScores;
                    double ionMass = peak.Mass;
                    double fragmentMass = ionMass - ionType.OffsetComposition.Mass;
                    var binMass = ionType.IsPrefix ? fragmentMass : proteinMass - fragmentMass;
                    var binIndex = this.massBins.GetBinNumber(binMass);
                    if (binIndex >= 0 && binIndex < numNodes)
                    {
                        var score = scorer.GetFragmentScore(peak, ionType); // Error score will come from the edge weight
                        if (scores[binIndex] == null || score >= scores[binIndex])
                        {
                            scores[binIndex] = score;
                        }
                    }
                }
            }

            return nodeScores;
        }
    }
}
