﻿using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring
{
    public class CorrMatchedPeakCounter : IScorer
    {
        public CorrMatchedPeakCounter(ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge, double corrScoreThreshold = 0.7)
        {
            _ms2Spec = ms2Spec;
            _tolerance = tolerance;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _corrScoreThreshold = corrScoreThreshold;
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

                var containsIon = false;
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var ion = new Ion(fragmentComposition, charge);
                    if (_ms2Spec.GetCorrScore(ion, _tolerance) > _corrScoreThreshold)
                    {
                        containsIon = true;
                        break;
                    }
                }

                if (containsIon) score += 1.0;
            }
            return score;
        }

        private readonly ProductSpectrum _ms2Spec;
        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly BaseIonType[] _baseIonTypes;
        private readonly double _corrScoreThreshold;

        public static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static CorrMatchedPeakCounter()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }
    }
}
