using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Scoring weights for entire feature vector, including all ions
    /// </summary>
    public class FeatureWeights
    {
        /// <summary>
        /// Gets or sets the intercept weight.
        /// </summary>
        public double Intercept { get; set; }

        /// <summary>
        /// Gets or sets the feature weights for each base ion type.
        /// </summary>
        public Dictionary<BaseIonType, IonWeights> IonWeights { get; set; }

        /// <summary>
        /// Gets or sets the weighting of mass error by the scoring model.
        /// </summary>
        public double MassErrorWeight { get; set; }

        /// <summary>
        /// Get the weighted score of the ion including mass error in the score.
        /// </summary>
        /// <param name="baseIonType">The base ion type to use the feature vector for.</param>
        /// <param name="peak">The peak for the ion to calculate the score for.</param>
        /// <param name="theoreticalMass">Theoretical mass</param>
        /// <returns>The score.</returns>
        public double GetMatchedIonPeakScoreWithError(BaseIonType baseIonType, DeconvolutedPeak peak, double theoreticalMass)
        {
            double score = featureVector.GetMatchedIonPeakScoreWithoutError(peak);
            var featureVector = IonWeights[baseIonType];

            var ppmError = Math.Abs(peak.Mass - theoreticalMass) / theoreticalMass * 1e6;
            score += MassErrorWeight * ppmError;

            return score;
        }
    }
}
