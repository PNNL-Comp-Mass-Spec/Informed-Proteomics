using System;
using System.Collections.Generic;

using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Factory for selecting the correct scoring parameters for a certain spectrum, activation method, and mass.
    /// </summary>
    public class FlipScorerFactory : IFragmentScorerFactory
    {
        /// <summary>
        /// For selecting deconvoluted spectra required for scoring.
        /// </summary>
        private readonly DPbfLcMsRun lcmsRun;

        /// <summary>
        /// For selecting the correct scoring parameters.
        /// </summary>
        private readonly ScoringParameterSet scoringParameterSet;

        /// <summary>
        /// Cached scorers.
        /// </summary>
        private readonly Dictionary<int, IScorer> scorerHash;

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScorerFactory" /> class with a preloaded scoring parameter set.
        /// </summary>
        /// <param name="lcmsRun">LcmsRun for selected deconvoluted spectra required for scoring.</param>
        /// <param name="scoringParameterSet">Existing scoring parameters to select from.</param>
        public FlipScorerFactory(DPbfLcMsRun lcmsRun, ScoringParameterSet scoringParameterSet)
        {
            this.lcmsRun = lcmsRun;
            this.scoringParameterSet = scoringParameterSet;
            this.scorerHash = new Dictionary<int, IScorer>();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScorerFactory" /> class with a new scoring parameter set.
        /// </summary>
        /// <param name="lcmsRun">LcmsRun for selected deconvoluted spectra required for scoring.</param>
        /// <param name="isTopDown">A value indicating whether the experiment is an intact experiment or not.</param>
        public FlipScorerFactory(DPbfLcMsRun lcmsRun, bool isTopDown = true) 
            : this(lcmsRun, new ScoringParameterSet(isTopDown)) { }

        /// <summary>
        /// Get the scorer by spectrum.
        /// This creates the scorer directly for that spectrum. It bypasses the scorer cache.
        /// </summary>
        /// <param name="spectrum">The spectrum number to get the scorer for.</param>
        /// <param name="precursorMass">The precursor mass to get the scorer for.</param>
        /// <param name="precursorCharge">The precursor charge to get the scorer for.</param>
        /// <returns>The scorer selected based on the arguments.</returns>
        public IScorer GetScorer(ProductSpectrum spectrum, double precursorMass, int precursorCharge)
        {
            var parameters = this.scoringParameterSet.GetScoringParameters(spectrum.ActivationMethod, precursorMass);
            var targetedDeconvolutedSpectrum = new TargetedDeconvolutedSpectrum(spectrum, 1, precursorCharge - 1);
            return new FlipScorer<TargetedDeconvolutedSpectrum>(parameters, targetedDeconvolutedSpectrum);
        }

        /// <summary>
        /// Get the scorer by scan number.
        /// </summary>
        /// <param name="scanNum">The scan number to get the scorer for.</param>
        /// <param name="precursorMass">The precursor mass to get the scorer for.</param>
        /// <param name="precursorCharge">The precursor charge to get the scorer for.</param>
        /// <param name="activationMethod">The activation method to get the scorer for.</param>
        /// <returns>The scorer selected based on the arguments.</returns>
        public IScorer GetScorer(int scanNum, double precursorMass = 0, int precursorCharge = 1, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            IScorer scorer;
            if (this.scorerHash.ContainsKey(scanNum))
            {
                scorer = this.scorerHash[scanNum];
            }
            else
            {
                var parameters = this.scoringParameterSet.GetScoringParameters(activationMethod, precursorMass);
                scorer = new FlipScorer<DeconvolutedSpectrum>(parameters, this.lcmsRun.GetSpectrum(scanNum) as DeconvolutedSpectrum);
                this.scorerHash.Add(scanNum, scorer);
            }

            return scorer;
        }
    }
}
