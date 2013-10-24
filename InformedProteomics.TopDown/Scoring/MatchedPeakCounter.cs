using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class MatchedPeakCounter: IScorer
    {
        public MatchedPeakCounter(Spectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge)
        {
            _ms2Spec = ms2Spec;
            _tolerance = tolerance;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
        }

        public double GetPrecursorIonScore(Ion precursorIon)
        {
            return 0;
        }

        public double GetFragmentScore(Ion precursorIon, Composition suffixFragmentComposition)
        {
            var score = 0.0;

            var prefixFragmentComposition = precursorIon.Composition - suffixFragmentComposition;
            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;
                fragmentComposition.ComputeApproximateIsotopomerEnvelop();

                if (fragmentComposition.GetMass() < 0)
                {
                    Console.WriteLine("************* {0}, {1}, {2}", suffixFragmentComposition, fragmentComposition, baseIonType.Symbol);
                }
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var ion = new Ion(fragmentComposition, charge);
                    if (_ms2Spec.ContainsIon(ion, _tolerance, RelativeIsotopeIntensityThreshold))
                    {
                        score += 1.0;
                    }
                    //var isotopes = ion.GetIsotopes(RelativeIsotopeIntensityThreshold);
                    //var allIonsExist = true;
                    //foreach (var isotope in isotopes)
                    //{
                    //    var isotopeIndex = isotope.Item1;
                    //    var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                    //    if (_ms2Spec.FindPeak(isotopeMz, _tolerance) == null)
                    //    {
                    //        allIonsExist = false;
                    //        break;
                    //    }
                    //}
                    //if (allIonsExist) score += 1.0;
                }
            }
            return score;
        }

        private readonly Spectrum _ms2Spec;
        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;

        private const double RelativeIsotopeIntensityThreshold = 0.7;

        private static readonly BaseIonType[] BaseIonTypes;
        static MatchedPeakCounter()
        {
            BaseIonTypes = new[] { BaseIonType.C, BaseIonType.Z };
        }

    }
}
