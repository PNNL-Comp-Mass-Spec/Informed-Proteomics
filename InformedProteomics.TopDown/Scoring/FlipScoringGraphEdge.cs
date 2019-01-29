using System;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.TopDown.Scoring.FlipScoring;

namespace InformedProteomics.TopDown.Scoring
{
    /// <summary>
    /// Scoring graph edge for a scoring graph that using the FLIP scoring model.
    /// </summary>
    public class FlipScoringGraphEdge : IScoringGraphEdge
    {
        /// <summary>
        /// The FLIP scoring model.
        /// </summary>
        private readonly FlipScorer<DeconvolutedSpectrum> scorer;

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScoringGraphEdge" />
        /// with initialized data.
        /// </summary>
        /// <param name="prevNodeIndex">The source node index.</param>
        /// <param name="sinkNodeIndex">The sink node index.</param>
        /// <param name="weight">The edge weight.</param>
        /// <param name="label">The amino acid that that the edge corresponds to.</param>
        /// <param name="scorer">The FLIP scoring model.</param>
        public FlipScoringGraphEdge(
                    int prevNodeIndex,
                    int sinkNodeIndex,
                    double weight,
                    AminoAcid label,
                    FlipScorer<DeconvolutedSpectrum> scorer)
        {
            PrevNodeIndex = prevNodeIndex;
            SinkNodeIndex = sinkNodeIndex;
            Weight = weight;
            Label = label;
            this.scorer = scorer;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScoringGraphEdge" />
        /// with uninitialized data.
        /// </summary>
        /// <param name="scorer">The FLIP scoring model.</param>
        public FlipScoringGraphEdge(FlipScorer<DeconvolutedSpectrum> scorer)
            : this(0, 0, 0, null, scorer) { }

        /// <summary>
        /// Gets or sets the source node index.
        /// </summary>
        public int PrevNodeIndex { get; set; }

        /// <summary>
        /// Gets or sets the sink node index.
        /// </summary>
        public int SinkNodeIndex { get; set; }

        /// <summary>
        /// Gets or sets the weight of the node.
        /// </summary>
        public double Weight { get; set; }

        /// <summary>
        /// Gets or sets the amino acid that the edge corresponds to.
        /// </summary>
        public AminoAcid Label { get; set; }

        /// <summary>
        /// Gets of the score of the edge.
        /// </summary>
        /// <returns>The score of the edge.</returns>
        public double GetScore()
        {
            var massError = Math.Abs(SinkNodeIndex - (PrevNodeIndex + Label.Mass)) / SinkNodeIndex * 1e6;
            var errorScore = scorer.GetErrorScore(massError);

            return errorScore + Weight;
        }
    }
}
