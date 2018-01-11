using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MathAndStats;

namespace InformedProteomics.Tests.DevTests.TopDownAnalysis
{
    public class AnalysisTopDownMatchedPeaks
    {
        public const int MinCharge = 1;
        public const int MaxCharge = 20;

        public void OutputStatistics(ProductSpectrum spectrum, Sequence sequence)
        {
            var baseIonTypes = spectrum.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCid : BaseIonTypesEtd;
            var cleavages = sequence.GetInternalCleavages().ToArray();
            var tolerance = new Tolerance(10);

            var maxIntensity = spectrum.Peaks.Max(p => p.Intensity);

            foreach (var c in cleavages)
            {
                foreach (var baseIonType in baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                        ? c.PrefixComposition + baseIonType.OffsetComposition
                        : c.SuffixComposition + baseIonType.OffsetComposition;

                    for (var charge = MinCharge; charge <= MaxCharge; charge++)
                    {
                        var ion = new Ion(fragmentComposition, charge);
                        var observedPeaks = spectrum.GetAllIsotopePeaks(ion, tolerance, RelativeIsotopeIntensityThreshold);

                        if (observedPeaks == null) continue;

                        var mostAbundantIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();

                        // representative peak intensity
                        var ionPeakIntensity = observedPeaks[mostAbundantIsotopeIndex].Intensity;

                        // calc. correlation
                        var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
                        var observedIntensities = new double[observedPeaks.Length];
                        for (var i = 0; i < observedPeaks.Length; i++)
                        {
                            var observedPeak = observedPeaks[i];
                            observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                        }
                        var corrCoeff = FitScoreCalculator.GetPearsonCorrelation(isotopomerEnvelope, observedIntensities);

                        // mz error
                        var mostAbundantIsotopeMz = ion.GetIsotopeMz(mostAbundantIsotopeIndex);
                        var errorPpm = ((observedPeaks[mostAbundantIsotopeIndex].Mz - mostAbundantIsotopeMz)/
                                        mostAbundantIsotopeMz)*1e6;
                    }
                }
            }
        }

        private const double RelativeIsotopeIntensityThreshold = 0.7;
        public static readonly BaseIonType[] BaseIonTypesCid, BaseIonTypesEtd;

        static AnalysisTopDownMatchedPeaks()
        {
            BaseIonTypesCid = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesEtd = new[] { BaseIonType.C, BaseIonType.Z };
        }
    }
}
