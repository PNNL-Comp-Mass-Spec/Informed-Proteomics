using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.Interfaces;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Top down scorer that uses the FLIP scoring model.
    /// </summary>
    public class FlipScorer<TNeutralMonoSpectrum> : IInformedScorer, IScorer
        where TNeutralMonoSpectrum : IDeconvolutedSpectrum
    {
        /// <summary>
        /// Parameters to use for scoring ions.
        /// </summary>
        private readonly ScoringParameters scoringParameters;

        /// <summary>
        /// Initializes a new instance of the <see cref="FlipScorer{TNeutralMonoSpectrum}" /> class.
        /// </summary>
        /// <param name="scoringParameters">Parameters to use for scoring ions.</param>
        /// <param name="productSpectrum">The deconvoluted spectrum to get ions to calculate scores for.</param>
        public FlipScorer(ScoringParameters scoringParameters, TNeutralMonoSpectrum productSpectrum)
        {
            this.scoringParameters = scoringParameters;
            this.ProductSpectrum = productSpectrum;
        }

        /// <summary>
        /// Gets the the ions selected for scoring for this bin.
        /// </summary>
        public BaseIonType[] SelectedIonTypes => this.scoringParameters.SelectedIonTypes;

        /// <summary>
        /// Gets or sets the deconvoluted spectrum to get ions to calculate scores for.
        /// </summary>
        public TNeutralMonoSpectrum ProductSpectrum { get; private set; }

        /// <summary>
        /// Gets the score at sequence cleavage summed for all ions that were selected during training.
        /// </summary>
        /// <param name="nTerminalFragmentComposition">Composition of NTerminal sequence cleavage.</param>
        /// <param name="cTerminalFragmentComposition">Composition of CTerminal sequence cleavage.</param>
        /// <param name="nTerminalResidue">Residue at NTerminal sequence cleavage.</param>
        /// <param name="cTerminalResidue">Residue at CTerminal sequence cleavage.</param>
        /// <returns>The summed score of all ions.</returns>
        public double GetFragmentScore(
            Composition nTerminalFragmentComposition,
            Composition cTerminalFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
        {
            double score = this.scoringParameters.FeatureWeights.Intercept;
            foreach (BaseIonType ionType in this.scoringParameters.SelectedIonTypes)
            {
                FeatureWeights featureVector = this.scoringParameters.FeatureWeights;
                Composition composition = ionType.IsPrefix ? nTerminalFragmentComposition : cTerminalFragmentComposition;
                AminoAcid residue = ionType.IsPrefix ? nTerminalResidue : cTerminalResidue;
                IEnumerable<Composition> offsetCompositions = ionType.GetPossibleCompositions(residue);

                double ionScore = 0.0;
                foreach (Composition offsetComp in offsetCompositions)  // Iterate over possible compositions for the current base ion type
                {
                    var ionComp = composition + offsetComp;
                    var peak = this.ProductSpectrum.FindPeak(ionComp, this.scoringParameters.TrainerTolerance);
                    if (peak == null)
                    {
                        continue;
                    }

                    // Use only the score of the best scoring composition option (we have all the beset scoring compositions!)
                    var compScore = featureVector.GetMatchedIonPeakScoreWithError(ionType, peak, ionComp.Mass);
                    if (compScore > ionScore)
                    {
                        ionScore = compScore;
                    }
                }

                score += ionScore;
            }

            return score;
        }

        /// <summary>
        /// Gets the score for a single deconvoluted peak.
        /// Does not include Mass error score.
        /// </summary>
        /// <param name="deconvolutedPeak">The peak to score.</param>
        /// <param name="ionType">The type of ion to score against.</param>
        /// <returns>The score of the peak.</returns>
        public double GetFragmentScore(DeconvolutedPeak deconvolutedPeak, BaseIonType ionType)
        {
            IonWeights featureVector = this.scoringParameters.FeatureWeights.IonWeights[ionType];
            return featureVector.GetMatchedIonPeakScoreWithoutError(deconvolutedPeak);
        }

        /// <summary>
        /// Gets the mass error score for a given mass error.
        /// </summary>
        /// <param name="massError">The mass error to use for scoring.</param>
        /// <returns>The </returns>
        public double GetErrorScore(double massError)
        {
            return massError * this.scoringParameters.FeatureWeights.MassErrorWeight;
        }

        /// <summary>
        /// Gets the number of sequence fragments found in the spectrum.
        /// </summary>
        /// <param name="sequence">The sequence to count fragments for.</param>
        /// <returns>The number of matched fragments.</returns>
        public int GetNumMatchedFragments(Sequence sequence)
        {
            int count = 0;
            var cleavages = sequence.GetInternalCleavages();
            foreach (var cl in cleavages)
            {
                var fragScore = this.GetFragmentScore(cl.PrefixComposition, cl.SuffixComposition, cl.PrefixResidue, cl.SuffixResidue);
                count += !fragScore.Equals(0.0) ? 1 : 0;
            }

            return count;
        }

        /// <summary>
        /// Gets the user-presentable score that will be displayed on results file.
        /// This is the probability based on the logistic regression scoring.
        /// </summary>
        /// <param name="sequence">The sequence to calculate a score for.</param>
        /// <returns>The score.</returns>
        public double GetUserVisibleScore(Sequence sequence)
        {
            // Calculate linear combination.
            double score = 0.0;
            var cleavages = sequence.GetInternalCleavages();
            foreach (var cl in cleavages)
            {
                score += this.GetFragmentScore(cl.PrefixComposition, cl.SuffixComposition, cl.PrefixResidue, cl.SuffixResidue);
            }

            // Calculate probability
            var eta = Math.Exp(this.scoringParameters.FeatureWeights.Intercept + score);
            return eta / (eta + 1);
        }

        /// <summary>
        /// Gets the score cut off for filtering out PSMs before they are passed into the generating function.
        /// </summary>
        public double ScoreCutOff { get { return 0.0; } }
    }
}
