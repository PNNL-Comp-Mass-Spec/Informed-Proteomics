using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Deconvoluter
    {
        public static DeconvolutedSpectrum GetDeconvolutedSpectrum(
                    Spectrum spec, int minCharge, int maxCharge,
                    int isotopeOffsetTolerance, double filteringWindowSize,
                    Tolerance tolerance, double corrScoreThreshold = 0.7)
        {
            var peaks = GetDeconvolutedPeaks(spec.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize,
                tolerance, corrScoreThreshold);

            return new DeconvolutedSpectrum(spec, peaks.ToArray());
        }

        public static DeconvolutedSpectrum GetCombinedDeconvolutedSpectrum(
            Spectrum spec,
            int minCharge,
            int maxCharge,
            int isotopeOffsetTolerance,
            double filteringWindowSize,
            Tolerance tolerance,
            double corrScoreThreshold)
        {
            var deconvolutedPeaks = GetDeconvolutedPeaks(spec, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize, tolerance, corrScoreThreshold);
            var binHash = new Dictionary<int, DeconvolutedPeak>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = (int)Math.Round(mass * Constants.RescalingConstantHighPrecision);
                if (!binHash.ContainsKey(binNum))
                {
                    binHash.Add(binNum, deconvolutedPeak);
                }
                else if (binHash[binNum].Intensity < deconvolutedPeak.Intensity)
                {
                    binHash[binNum] = deconvolutedPeak;
                }
            }

            var peakList = binHash.Values.OrderBy(peak => peak.Mz).ToArray();
            return new DeconvolutedSpectrum(spec, peakList);
        }

        public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
            Spectrum spec, int minCharge, int maxCharge,
            int isotopeOffsetTolerance, double filteringWindowSize,
            Tolerance tolerance, double corrScoreThreshold)
        {
            return GetDeconvolutedPeaks(spec.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize,
                tolerance, corrScoreThreshold);
        }

        // Select the best peak within +/- filteringWindowSize
        public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
            Peak[] peaks, int minCharge, int maxCharge, 
            int isotopeOffsetTolerance, double filteringWindowSize,
            Tolerance tolerance, double corrScoreThreshold)
        {

            var monoIsotopePeakList = new List<DeconvolutedPeak>();
            for (var peakIndex = 0; peakIndex < peaks.Length; peakIndex++)
            {
                var peak = peaks[peakIndex];

                // Check whether peak has the maximum intensity within the window
                var isBest = true;

                var prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if ((peak.Mz - prevPeak.Mz) > filteringWindowSize) break;
                    if (prevPeak.Intensity > peak.Intensity)
                    {
                        isBest = false;
                        break;
                    }
                    prevIndex--;
                }

                if (!isBest) continue;

                var nextIndex = peakIndex + 1;
                while (nextIndex < peaks.Length)
                {
                    var nextPeak = peaks[nextIndex];
                    if ((nextPeak.Mz - peak.Mz) > filteringWindowSize) break;
                    if (nextPeak.Intensity > peak.Intensity)
                    {
                        isBest = false;
                        break;
                    }
                    nextIndex++;
                }

                if (!isBest) continue;

                // peak has the maximum intensity, window = [prevIndex+1,nextIndex-1]

                var window = new Peak[nextIndex - prevIndex - 1];
                Array.Copy(peaks, prevIndex + 1, window, 0, window.Length);
                var windowSpectrum = new Spectrum(window, 1);
                var peakMz = peak.Mz;

                double bestCorrelation = 0.0;
                DeconvolutedPeak bestPeak = null;

                for (var charge = maxCharge; charge >= minCharge; charge--)
                {
                    var mass = peak.Mz * charge;
                    var mostAbundantIsotopeIndex = Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex;

                    for (var isotopeIndex = mostAbundantIsotopeIndex - isotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex + isotopeOffsetTolerance; isotopeIndex++)
                    {
                        var monoIsotopeMass = Ion.GetMonoIsotopicMass(peakMz, charge, isotopeIndex);
                        var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(monoIsotopeMass);
                        var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
                        if (observedPeaks == null) continue;

                        var envelop = isotopomerEnvelope.Envolope;
                        var observedIntensities = new double[observedPeaks.Length];

                        for (var i = 0; i < observedPeaks.Length; i++)
                        {
                            var observedPeak = observedPeaks[i];
                            observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                        }

                        var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelop, observedIntensities);
                        var bcDist = sim.Item1;
                        var corr = sim.Item2;

                        if (corr < corrScoreThreshold && bcDist > 0.03) continue;

                        // monoIsotopeMass is valid
                        if (corr >= bestCorrelation)
                        {
                            bestCorrelation = corr;
                            bestPeak = new DeconvolutedPeak(monoIsotopeMass, observedIntensities[mostAbundantIsotopeIndex], charge, corr, bcDist, observedPeaks);
                        }
                    }
                }

                if (bestPeak != null)
                {
                    monoIsotopePeakList.Add(bestPeak);
                }
            }

            monoIsotopePeakList.Sort();
            return monoIsotopePeakList;
        }

        public static List<DeconvolutedPeak> FilterOut(List<DeconvolutedPeak> peaks, double windowSize, int topRankCutoff)
        {
            var retSpec = new List<DeconvolutedPeak>();
            for (var peakIndex = 0; peakIndex < peaks.Count; peakIndex++)
            {
                var rank = 1;
                var thisPeak = peaks[peakIndex];
                var thisMass = thisPeak.Mass;
                var thisInten = thisPeak.Intensity;

                // move left
                var prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if (thisMass - prevPeak.Mass > windowSize) break;
                    if (prevPeak.Intensity > thisInten) rank++;
                    prevIndex--;
                }

                // move right
                var nextIndex = peakIndex + 1;
                while (nextIndex < peaks.Count)
                {
                    var nextPeak = peaks[nextIndex];
                    if (nextPeak.Mass - thisMass > windowSize) break;
                    if (nextPeak.Intensity > thisInten) rank++;
                    nextIndex++;
                }
                if (rank <= topRankCutoff) retSpec.Add(thisPeak);
            }                

            return retSpec;
        }

        
    }
}
