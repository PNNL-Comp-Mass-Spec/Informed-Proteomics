using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class LikelihoodScorer : IScorer
    {
        public const double RelativeIntensityThreshold = 0.1;
        public LikelihoodScorer(LikelihoodScoringModel model, ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge)
        {
            _model = model;
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

                var bestCosineScore = 0.0;
                var bestObsIntensity = -1.0;
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var ion = new Ion(fragmentComposition, charge);

                    var observedPeaks = _ms2Spec.GetAllIsotopePeaks(ion, _tolerance, RelativeIntensityThreshold);
                    if (observedPeaks == null) continue;

                    var theoIntensities = new double[observedPeaks.Length];
                    var observedIntensities = new double[observedPeaks.Length];

                    for (var i = 0; i < observedPeaks.Length; i++)
                    {
                        theoIntensities[i] = isotopomerEnvelope[i];
                        var observedPeak = observedPeaks[i];
                        observedIntensities[i] = observedPeak != null ? observedPeak.Intensity : 0.0;
                    }
                    var cosineScore = FitScoreCalculator.GetCosine(isotopomerEnvelope, observedIntensities);
                    if (cosineScore > bestCosineScore)
                    {
                        bestCosineScore = cosineScore;
                        bestObsIntensity = observedIntensities.Max();
                    }
                }
//                Console.WriteLine("{0} {1}: {2}", baseIonType.Symbol, prefixFragmentComposition, _model.GetScore(baseIonType, bestCosineScore, bestObsIntensity));
                score += _model.GetScore(baseIonType, bestCosineScore, bestObsIntensity);
            }
            return score;
        }

        private readonly LikelihoodScoringModel _model;
        private readonly ProductSpectrum _ms2Spec;
        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly BaseIonType[] _baseIonTypes;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static LikelihoodScorer()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

    }
}
