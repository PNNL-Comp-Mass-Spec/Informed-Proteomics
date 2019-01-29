using System;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.MathAndStats;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// This class exists to do a single targeted, single ion deconvolution of a non-deconvoluted spectrum.
    /// </summary>
    public class TargetedDeconvolutedSpectrum : IDeconvolutedSpectrum
    {
        /// <summary>
        /// Non-deconvoluted spectrum.
        /// </summary>
        private readonly Spectrum spectrum;

        /// <summary>
        /// The minimum charge state to consider.
        /// </summary>
        private readonly int minCharge;

        /// <summary>
        /// The maximum charge state to consider.
        /// </summary>
        private readonly int maxCharge;

        /// <summary>
        /// Initializes a new instance of the <see cref="TargetedDeconvolutedSpectrum" /> class.
        /// </summary>
        /// <param name="spectrum">Non-deconvoluted spectrum.</param>
        /// <param name="minCharge">The minimum charge state to consider.</param>
        /// <param name="maxCharge">The maximum charge state to consider.</param>
        public TargetedDeconvolutedSpectrum(Spectrum spectrum, int minCharge = 1, int maxCharge = 20)
        {
            this.spectrum = spectrum;
            this.minCharge = minCharge;
            this.maxCharge = maxCharge;
        }

        /// <summary>
        /// Find a neutral monoisotopic peak corresponding to a single compound
        /// </summary>
        /// <param name="composition">The compound to look for.</param>
        /// <param name="tolerance">The peak tolerance to use.</param>
        /// <returns>The neutral monoisotopic peak.</returns>
        public DeconvolutedPeak FindPeak(Composition.Composition composition, Tolerance tolerance)
        {
            var mass = composition.Mass;
            DeconvolutedPeak deconvolutedPeak = null;
            var maxIntensity = 0.0;
            for (var charge = minCharge; charge <= maxCharge; charge++)
            {
                var ion = new Ion(composition, charge);
                var isotopePeaks = spectrum.GetAllIsotopePeaks(ion, tolerance);
                var intensity = isotopePeaks.Max(peak => peak.Intensity);
                var corrCos = GetCorrCos(ion, isotopePeaks);
                if (intensity > maxIntensity)
                {
                    deconvolutedPeak = new DeconvolutedPeak(mass, intensity, charge, corrCos.Item1, corrCos.Item2, isotopePeaks);
                }
            }

            return deconvolutedPeak;
        }

        /// <summary>
        /// Gets the pearson correlation and cosine score of an observed isotopic distribution compared to a theoretical ion.
        /// </summary>
        /// <param name="ion">Theoretical ion.</param>
        /// <param name="observedPeaks">The observed isotopic distribution.</param>
        /// <returns>A tuple where the first item is pearson correlation and the second item is cosine.</returns>
        private Tuple<double, double> GetCorrCos(Ion ion, Peak[] observedPeaks)
        {
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
            var envelope = isotopomerEnvelope;
            var observedIntensities = new double[envelope.Length];
            for (var i = 0; i < isotopomerEnvelope.Length; i++)
            {
                if (observedPeaks[i] != null) observedIntensities[i] = observedPeaks[i].Intensity;
            }

            var pearsonCorrelation = FitScoreCalculator.GetPearsonCorrelation(envelope, observedIntensities);
            var cosine = FitScoreCalculator.GetCosine(envelope, observedIntensities);

            return new Tuple<double, double>(pearsonCorrelation, cosine);
        }
    }
}
