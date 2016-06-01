using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring
{
    public class MatchedPeakCounter: AbstractFragmentScorer
    {
        public MatchedPeakCounter(ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge)
            : base(ms2Spec, tolerance, minCharge, maxCharge)
        {
            RelativeIsotopeIntensityThreshold = 0.7;
        }

        public override double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
        {
            var score = 0.0;

            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;

                if (fragmentComposition.Mass < Ms2Spectrum.Peaks[0].Mz) continue;
                var chargeRange = GetMinMaxChargeRange(fragmentComposition);

                var containsIon = false;
                for (var charge = chargeRange.Min; charge <= chargeRange.Max; charge++)
                {
                    var ion = new Ion(fragmentComposition, charge);
                    if (Ms2Spectrum.ContainsIon(ion, Tolerance, RelativeIsotopeIntensityThreshold))
                    {
                        containsIon = true;
                        break;
                    }
                }

                if (containsIon) score += 1.0;
            }
            return score;
        }
    }
}
