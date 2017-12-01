using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// The score weights for a single base ion type.
    /// </summary>
    public class IonWeights
    {
        /// <summary>
        /// Gets or sets the weighting given to the summed intensity of the ion.
        /// </summary>
        public double SummedIntensityWeight { get; set; }

        /// <summary>
        /// Gets or sets the weighting given to the summed pearson correlation of the ion.
        /// </summary>
        public double SummedPearsonCorrelationWeight { get; set; }

        /// <summary>
        /// Gets or sets the weighting given to the summed cosine of the ion.
        /// </summary>
        public double SummedCosineWeight { get; set; }

        /// <summary>
        /// Get the weighted score of the ion without including mass error in the score.
        /// </summary>
        /// <param name="peak">The peak for the ion to calculate the score for.</param>
        /// <returns>The score.</returns>
        public double GetMatchedIonPeakScoreWithoutError(DeconvolutedPeak peak)
        {
            double score = 0.0;
            score += this.SummedIntensityWeight * peak.Intensity;
            score += this.SummedPearsonCorrelationWeight * peak.Corr;
            score += this.SummedCosineWeight * peak.Dist;
            return score;
        }
    }
}
