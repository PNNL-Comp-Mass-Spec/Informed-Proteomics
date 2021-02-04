using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.TopDown.Scoring.FlipScoring;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    /// <summary>
    /// Encapsulates the logic of creating the correct scoring graph for a certain scorer.
    /// </summary>
    public class ScoringGraphFactory
    {
        /// <summary>
        /// The mass bins for constructing the nodes of the graph.
        /// </summary>
        private readonly IMassBinning massBins;

        /// <summary>
        /// Amino acid set to build the graph edges from.
        /// </summary>
        private readonly AminoAcidSet aminoAcidSet;

        /// <summary>
        /// The amino acid probabilities.
        /// </summary>
        private readonly Dictionary<char, double> aminoAcidProbabilities;

        /// <summary>
        /// Initializes a new instance of the <see cref="ScoringGraphFactory" /> class.
        /// </summary>
        /// <param name="massBins">The mass bins for constructing the nodes of the graph.</param>
        /// <param name="aminoAcidSet">Amino acid set to build the graph edges from.</param>
        /// <param name="aminoAcidProbabilities">The amino acid probabilities.</param>
        public ScoringGraphFactory(IMassBinning massBins, AminoAcidSet aminoAcidSet, Dictionary<char, double> aminoAcidProbabilities = null)
        {
            this.massBins = massBins;
            this.aminoAcidSet = aminoAcidSet;
            this.aminoAcidProbabilities = aminoAcidProbabilities;
        }

        /// <summary>
        /// Create the correct scoring graph for the given scorer and parent mass.
        /// </summary>
        /// <param name="scorer">The scorer to create the scoring graph for.</param>
        /// <param name="parentMass">The maximum mass of the scoring graph.</param>
        /// <returns>The initialized scoring graph.</returns>
        public IScoringGraph GetScoringGraph(IScorer scorer, double parentMass)
        {
            if (scorer is CompositeScorerBasedOnDeconvolutedSpectrum compositeScorer)
            {
                var scoringGraphFactory = new ProteinScoringGraphFactory(this.massBins, this.aminoAcidSet);
                return scoringGraphFactory.CreateScoringGraph(compositeScorer, parentMass);
            }
            else if (scorer is FlipScorer<DeconvolutedSpectrum> flipScorer)
            {
                var scoringGraphFactory = new FlipScoringGraphFactory(this.massBins, this.aminoAcidSet, this.aminoAcidProbabilities);
                return scoringGraphFactory.GetScoringGraph(flipScorer, parentMass);
            }
            else
            {
                throw new ArgumentException("No scoring graph exists for that scorer.");
            }
        }
    }
}
