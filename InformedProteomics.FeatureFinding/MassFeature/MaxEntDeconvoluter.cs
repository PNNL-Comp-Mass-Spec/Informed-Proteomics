using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public class MaxEntDeconvoluter
    {
        /// <summary>
        /// MaxEnt deconvolution algorithm constructor
        /// </summary>
        /// <param name="tolerance">tolerance</param>
        /// <param name="massBinning">mass binning interface</param>
        /// <param name="minCharge">maximum charge to be considered</param>
        /// <param name="maxCharge">minimum charge to be considered</param>
        /// <param name="minMass">minimum mass to be considered</param>
        /// <param name="maxMass">maximum mass to be considered</param>
        public MaxEntDeconvoluter(Spectrum spec, Tolerance tolerance, MzComparerWithBinning massBinning, int minCharge = 1, int maxCharge = 100, double minMass = 10000, double maxMass = 100000)
        {
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _minMass = minMass;
            _maxMass = maxMass;
            _tolerance = tolerance;
            _massBinning = massBinning;
            _spectrum = spec;
            _deconvolutedSpectrum = null;
        }

        /// <summary>
        /// Deconvolute mass spectrum using deconvolution algorithm described in Mann et al., Anal. Chem. 1989
        /// </summary>
        /// <param name="spec">spectrum</param>
        /// <returns>deconvoluted spectrum</returns>
        public DeconvolutedSpectrum GetDeconvolutedSpectrum()
        {
            if (_deconvolutedSpectrum != null) return _deconvolutedSpectrum;

            var minBinIndex = _massBinning.GetBinNumber(_minMass);
            var maxBinIndex = _massBinning.GetBinNumber(_maxMass);
            var massIntensity = new double[maxBinIndex - minBinIndex + 1];

            foreach (var peak in _spectrum.Peaks)
            {
                for (var charge = _minCharge; charge <= _maxCharge; charge++)
                {
                    var mass = (peak.Mz - Constants.Proton) * charge;
                    if (mass < _minMass) continue;
                    if (mass > _maxMass) break;

                    var massBinNum = _massBinning.GetBinNumber(mass);
                    massIntensity[massBinNum - minBinIndex] += peak.Intensity;
                }
            }
            var deconvPeaks = new List<DeconvolutedPeak>();
            for (var i = 0; i < massIntensity.Length; i++)
            {
                var intensity = massIntensity[i];
                if (intensity < float.Epsilon) continue;
                var mass = _massBinning.GetMzAverage(i + minBinIndex);
                deconvPeaks.Add(new DeconvolutedPeak(mass, intensity, 1, 0, 0));
            }
            _deconvolutedSpectrum = new DeconvolutedSpectrum(_spectrum, deconvPeaks.ToArray());
            return _deconvolutedSpectrum;
        }

        /// <summary>
        /// Finds mass features from deconvoluted spectrum. Monoiosotopic mass is estimated using averaging algorithm described in Mann et al., Anal. Chem. 1989
        /// </summary>
        /// <param name="deconvSpectrum">deconvoluted spectrum</param>
        /// <returns>mass features sorted by their abundancy</returns>
        public IEnumerable<Ms1Feature> GetMassFeatures()
        {
            var deconvSpectrum = GetDeconvolutedSpectrum();
            var deconvPeaks = (DeconvolutedPeak[])deconvSpectrum.Peaks;

            foreach (var deconvPeak in deconvPeaks.OrderByDescending(p => p.Intensity))
            {
                var collectedPeaks = CollectMultiplyChargedPeaks(deconvPeak.Mass, _spectrum);

                if (collectedPeaks.Count < 5) continue;

                var mass = GetAveragingMass(collectedPeaks);
                var charges = collectedPeaks.Select(p => p.Charge).ToArray();
                var abundance = collectedPeaks.Sum(p => p.Intensity);

                /*
                var maxChg = charges.Max();
                var minChg = charges.Min();
                var chargeDist = new double[maxChg - minChg + 1];
                foreach (var peak in collectedPeaks)
                {
                    chargeDist[peak.Charge - minChg] += peak.Intensity;
                }

                var mnStd = chargeDist.MeanStandardDeviation();
                var md = chargeDist.Median();
                var skewness = 3*(mnStd.Item1 - md)/md;
                if (skewness < -5 || skewness > 5) continue;
                */
                //yield return new Ms1Feature(_spectrum.ScanNum, mass, charges, abundance, collectedPeaks);
                yield return new Ms1Feature(_spectrum.ScanNum, mass, charges, abundance);
            }
        }

        /// <summary>
        /// Estimate monoisotopic mass from a set of peaks using averaging algorithm described in Mann et al., Anal. Chem. 1989
        /// </summary>
        /// <param name="peakSet">set of peaks</param>
        /// <returns>averaging monoisotopic mass</returns>
        public double GetAveragingMass(IList<DeconvolutedPeak> peakSet)
        {
            var weights = GetWeightingFactor(peakSet);
            var W = 0d;
            var mass = 0d;

            for (var i = 0; i < peakSet.Count; i++)
            {
                var iPeak = peakSet[i];
                var w = weights[i];
                mass += w * iPeak.Mass;
                W += w;
            }
            mass /= W;

            return mass;
        }

        private List<DeconvolutedPeak> CollectMultiplyChargedPeaks(double mass, Spectrum spectrum)
        {
            var deconvPeaks = new List<DeconvolutedPeak>();
            var minMz = spectrum.Peaks.First().Mz;
            var maxMz = spectrum.Peaks.Last().Mz;

            var minCharge = (int)Math.Max(_minCharge, Math.Ceiling(mass / (maxMz - Constants.Proton)));
            var maxCharge = (int)Math.Min(_maxCharge, Math.Floor(mass / (minMz - Constants.Proton)));

            for (var charge = minCharge; charge < maxCharge; charge++)
            {
                var mz = mass / charge + Constants.Proton;
                var peak = spectrum.FindPeak(mz, _tolerance);
                if (peak == null) continue;
                deconvPeaks.Add(new DeconvolutedPeak(peak, charge));
            }

            return deconvPeaks;
        }

        private double[] GetWeightingFactor(IList<DeconvolutedPeak> peakSet)
        {
            var maxCharge = peakSet.Max(p => p.Charge);
            var minCharge = peakSet.Min(p => p.Charge);
            var weights = new double[peakSet.Count];

            if (maxCharge == minCharge)
            {
                for (var i = 0; i < peakSet.Count; i++) weights[i] = 1;
                return weights;
            }

            for (var i = 0; i < peakSet.Count; i++)
            {
                var iPeak = peakSet[i];
                for (var j = 0; j < peakSet.Count; j++)
                {
                    var jPeak = peakSet[j];
                    if (iPeak.Charge == jPeak.Charge) continue;
                    var denominator = (double)iPeak.Charge / jPeak.Charge;
                    var numerator = Math.Abs(denominator - jPeak.MzWithoutAdductIonMass / iPeak.MzWithoutAdductIonMass);
                    weights[i] += Math.Pow(numerator / denominator, WeightingIndex);

                }
                weights[i] *= (maxCharge - minCharge);
            }

            return weights;
        }


        private const double WeightingIndex = 2;
        private readonly Tolerance _tolerance;
        private readonly double _minMass;
        private readonly double _maxMass;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly MzComparerWithBinning _massBinning;
        private readonly Spectrum _spectrum;
        private DeconvolutedSpectrum _deconvolutedSpectrum;

    }
}
