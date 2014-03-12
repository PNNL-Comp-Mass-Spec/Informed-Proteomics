using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class FitBasedLogLikelihoodRatioScorer : IScorer
    {
        public const double RelativeIntensityThreshold = 0.1;
       public FitBasedLogLikelihoodRatioScorer(ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge)
        {
            _ms2Spec = ms2Spec;
            _tolerance = tolerance;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _baseIonTypes = ms2Spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
        }

        public double GetPrecursorIonScore(Ion precursorIon)
        {
            return 0;
        }

        public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
        {
            var score = 0.0;

            foreach (var baseIonType in _baseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;
                fragmentComposition.ComputeApproximateIsotopomerEnvelop();
                var isotopomerEnvelope = fragmentComposition.GetIsotopomerEnvelop();

                //var bestFitScore = 1.0;
                var bestCosineScore = 0.0;
                var bestObsIntensity = -1.0;
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var ion = new Ion(fragmentComposition, charge);
                    var cosineScore = _ms2Spec.GetConsineScore(ion, _tolerance, RelativeIntensityThreshold);
                    if (cosineScore > bestCosineScore) bestCosineScore = cosineScore;
                    //var observedPeaks = _ms2Spec.GetAllIsotopePeaks(ion, _tolerance, RelativeIntensityThreshold);

                    //if (observedPeaks == null) continue;

                    //var theoIntensities = new float[observedPeaks.Length];
                    //var observedIntensities = new float[observedPeaks.Length];
                    //var maxObservedIntensity = float.NegativeInfinity;
                    //for (var i = 0; i < observedPeaks.Length; i++)
                    //{
                    //    theoIntensities[i] = isotopomerEnvelope[i];
                    //    if (observedPeaks[i] != null)
                    //    {
                    //        var observedIntensity = (float)observedPeaks[i].Intensity;
                    //        observedIntensities[i] = observedIntensity;
                    //        if (observedIntensity > maxObservedIntensity) maxObservedIntensity = observedIntensity;
                    //    }
                    //    else
                    //    {
                    //        observedIntensities[i] = 0;
                    //    }
                    //}

                    //for (var i = 0; i < observedPeaks.Length; i++)
                    //{
                    //    observedIntensities[i] /= maxObservedIntensity;
                    //}
                    //var fitScore = FitScoreCalculator.GetFitOfNormalizedVectors(isotopomerEnvelope, observedIntensities);
                    //if (fitScore < bestFitScore)
                    //{
                    //    bestFitScore = fitScore;
                    //    bestObsIntensity = maxObservedIntensity;
                    //}
                }
                //score += GetScore(baseIonType, bestFitScore, bestObsIntensity);
                score += GetScore(baseIonType, bestCosineScore, bestObsIntensity);
            }
            return score;
        }

        private readonly ProductSpectrum _ms2Spec;
        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly BaseIonType[] _baseIonTypes;

        private const double RelativeIsotopeIntensityThreshold = 0.8;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static FitBasedLogLikelihoodRatioScorer()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

        private static double GetScore(BaseIonType baseIonType, double fitScore, double intensity)
        {
            if (intensity < 0) return -1;
            return 0;
        }
    }
}
