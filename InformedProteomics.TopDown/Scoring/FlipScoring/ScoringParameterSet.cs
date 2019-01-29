using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Selects/reads the correct scoring parameters based on a scoring parameter description.
    /// </summary>
    public class ScoringParameterSet
    {
        /// <summary>
        /// A value indicating whether the experiment is an intact experiment or not.
        /// </summary>
        private readonly bool isTopDown;

        /// <summary>
        /// Cached scoring parameters.
        /// </summary>
        private readonly Dictionary<ScoringParameterDescription, ScoringParameters[]> scoringParameters;

        /// <summary>
        /// Initializes a new instance of the <see cref="ScoringParameterSet" /> class.
        /// </summary>
        /// <param name="isTopDown"></param>
        public ScoringParameterSet(bool isTopDown = true)
        {
            this.isTopDown = isTopDown;
            scoringParameters = new Dictionary<ScoringParameterDescription, ScoringParameters[]>();
        }

        /// <summary>
        /// Get the scoring parameter for the given activation method and mass bin.
        /// </summary>
        /// <param name="activationMethod">The activation method to find scoring parameters for.</param>
        /// <param name="precursorMass">The mass to find scoring param mass bin for.</param>
        /// <returns>The scoring parameters for the given activation method and mass bin.</returns>
        public ScoringParameters GetScoringParameters(ActivationMethod activationMethod, double precursorMass)
        {
            var parameters = GetScoringParameters(activationMethod);

            // Select bin edges for searching:
            var binEdges = parameters.Select(param => param.Mass).ToArray();

            // Find scoring parameters for precursor mass
            int binIndex = Array.BinarySearch(binEdges, precursorMass);
            if (binIndex < 0)
            {
                binIndex = ~binIndex;
                binIndex = Math.Min(binEdges.Length - 1, Math.Max(0, binIndex - 1));
            }

            return parameters[binIndex];
        }

        /// <summary>
        /// Gets the scoring parameter set for the given activation method.
        /// </summary>
        /// <param name="activationMethod">The activation method to find scoring parameters for.</param>
        /// <returns>Array of scoring parameters sorted by the mass of each parameter bin.</returns>
        public ScoringParameters[] GetScoringParameters(ActivationMethod activationMethod)
        {
            var searchParamDesc = new ScoringParameterDescription { ActivationMethod = activationMethod, IsTopDown = isTopDown };
            ScoringParameters[] parameters;
            if (scoringParameters.ContainsKey(searchParamDesc))
            {   // Seen this activation method already.
                parameters = scoringParameters[searchParamDesc];
            }
            else
            {   // Haven't seen this activation method yet. Try to load it from file.
                parameters = ScoringParameters.Parse(searchParamDesc.GetPath());
                scoringParameters.Add(searchParamDesc, parameters);
            }

            return parameters;
        }
    }
}
