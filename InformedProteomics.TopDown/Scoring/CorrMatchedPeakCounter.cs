using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring
{
    using InformedProteomics.Backend.Data.Sequence;

    public class CorrMatchedPeakCounter : AbstractFragmentScorer
    {
        public CorrMatchedPeakCounter(ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge, double corrScoreThreshold = 0.7):
            base(ms2Spec, tolerance, minCharge, maxCharge, 0.1)
        {
            _corrScoreThreshold = corrScoreThreshold;
        }

        public override double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
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
                    if (Ms2Spectrum.GetCorrScore(ion, Tolerance) > _corrScoreThreshold)
                    {
                        containsIon = true;
                        break;
                    }
                }

                if (containsIon) score += 1.0;
            }
            return score;
        }

        private readonly double _corrScoreThreshold;
    }
}
