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
            var peaks = GetDeconvolutedPeaks(spec.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, tolerance, corrScoreThreshold);

            return new DeconvolutedSpectrum(spec, peaks.ToArray());
        }

        //public static DeconvolutedSpectrum GetCombinedDeconvolutedSpectrum(
        //    Spectrum spec,
        //    int minCharge,
        //    int maxCharge,
        //    int isotopeOffsetTolerance,
        //    double filteringWindowSize,
        //    Tolerance tolerance,
        //    double corrScoreThreshold)
        //{
        //    var deconvolutedPeaks = GetDeconvolutedPeaks(spec, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize, tolerance, corrScoreThreshold);
        //    var binHash = new Dictionary<int, DeconvolutedPeak>();
        //    foreach (var deconvolutedPeak in deconvolutedPeaks)
        //    {
        //        var mass = deconvolutedPeak.Mass;
        //        var binNum = (int)Math.Round(mass * Constants.RescalingConstantHighPrecision);
        //        if (!binHash.ContainsKey(binNum))
        //        {
        //            binHash.Add(binNum, deconvolutedPeak);
        //        }
        //        else if (binHash[binNum].Intensity < deconvolutedPeak.Intensity)
        //        {
        //            binHash[binNum] = deconvolutedPeak;
        //        }
        //    }

        //    var peakList = binHash.Values.OrderBy(peak => peak.Mz).ToArray();
        //    return new DeconvolutedSpectrum(spec, peakList);
        //}

        public static DeconvolutedSpectrum GetCombinedDeconvolutedSpectrum(
            Spectrum spectrum,
            int minCharge,
            int maxCharge,
            int isotopeOffsetTolerance,
            Tolerance tolerance,
            double corrScoreThreshold)
        {
            var deconvolutedPeaks = GetDeconvolutedPeaks(spectrum.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, tolerance, corrScoreThreshold);
            var binHash = new Dictionary<int, DeconvolutedPeak>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = (int)Math.Round(mass * (2+Constants.RescalingConstant));
                if (!binHash.ContainsKey(binNum))
                {
                    binHash.Add(binNum, deconvolutedPeak);
                }
                else
                {
                    var binPeak = binHash[binNum];
                    if (binPeak.Intensity < deconvolutedPeak.Intensity)
                    {
                        binHash[binNum] = deconvolutedPeak;
                    }
                }
            }

            var peakList = binHash.Values.OrderBy(peak => peak.Mz).ToArray();
            return new DeconvolutedSpectrum(spectrum, peakList);
        }

        public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
            Peak[] peaks,
            int minCharge,
            int maxCharge,
            int isotopeOffsetTolerance,
            Tolerance tolerance,
            double corrScoreThreshold)
        {
            var spectrum = new Spectrum(peaks, 0);
            var monoIsotopePeakList = new List<DeconvolutedPeak>();

            var sortedPeaks = peaks.OrderByDescending(peak => peak.Intensity).ToArray();
            var peakUsed = new bool[peaks.Length];

            foreach (var peak in sortedPeaks)
            {
                var peakIndex = Array.BinarySearch(peaks, peak);
                if (peakUsed[peakIndex])
                {
                    continue;
                }

                double bestScore = 0.0;
                DeconvolutedPeak bestPeak = null;
                Tuple<Peak, int>[] bestObservedPeaks = null;

                for (int charge = minCharge; charge <= maxCharge; charge++)
                {
                    var mass = peak.Mz * charge - (charge * Constants.Proton);
                    var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(mass);
                    var mostAbundantIsotopeIndex = isotopomerEnvelope.MostAbundantIsotopeIndex;
                    var offsetTolerance = isotopeOffsetTolerance;
                    if (isotopeOffsetTolerance < 0)
                    {
                        offsetTolerance = isotopomerEnvelope.Envolope.Length;
                    }

                    for (var isotopeIndex = mostAbundantIsotopeIndex - offsetTolerance;
                         isotopeIndex <= mostAbundantIsotopeIndex + offsetTolerance;
                         isotopeIndex++)
                    {
                        var monoIsotopeMass = Ion.GetMonoIsotopicMass(peak.Mz, charge, isotopeIndex);

                        var observedPeaks = GetAllIsotopePeaks(spectrum, monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
                        if (observedPeaks == null) continue;

                        var envelop = isotopomerEnvelope.Envolope;
                        var observedIntensities = new double[observedPeaks.Length];

                        for (var i = 0; i < observedPeaks.Length; i++)
                        {
                            var observedPeak = observedPeaks[i];
                            if (observedPeak != null && peakUsed[observedPeak.Item2])
                            {
                                observedPeak = null;
                                observedPeaks[i] = null;
                            }

                            observedIntensities[i] = observedPeak != null ? (float)observedPeak.Item1.Intensity : 0.0;
                        }

                        var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelop, observedIntensities);
                        var bcDist = sim.Item1;
                        var corr = sim.Item2;

                        if (corr < corrScoreThreshold && bcDist > 0.03) continue;

                        var score = corr / (bcDist * (Math.Abs(mostAbundantIsotopeIndex - isotopeIndex) + 1));

                        //if (corr < corrScoreThreshold) continue;

                        // monoIsotopeMass is valid
                        if (score >= bestScore)
                        {
                            bestScore = score;
                            bestPeak = new DeconvolutedPeak(monoIsotopeMass, observedIntensities[mostAbundantIsotopeIndex], charge, corr, bcDist, observedPeaks.Where(p => p != null).Select(p => p.Item1).ToArray());
                            bestObservedPeaks = observedPeaks;
                        }
                    }
                }

                if (bestPeak != null)
                {
                    monoIsotopePeakList.Add(bestPeak);
                    foreach (var p in bestObservedPeaks)
                    {
                        if (p != null)
                        {
                            peakUsed[p.Item2] = true;   
                        }
                    }
                }
            }

            monoIsotopePeakList.Sort();
            return monoIsotopePeakList;
        }

        public static Tuple<Peak, int>[] GetAllIsotopePeaks(
            Spectrum spec,
            double monoisotopicMass,
            int charge,
            IsotopomerEnvelope envelope,
            Tolerance tolerance,
            double relativeIntensityThreshold)
        {
            var mostAbundantIsotopeIndex = envelope.MostAbundantIsotopeIndex;
            var isotopomerEnvelope = envelope.Envolope;
            var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, mostAbundantIsotopeIndex);
            var mostAbundantIsotopePeakIndex = spec.FindPeakIndex(mostAbundantIsotopeMz, tolerance);
            if (mostAbundantIsotopePeakIndex < 0) return null;

            var observedPeaks = new Tuple<Peak, int>[isotopomerEnvelope.Length];
            observedPeaks[mostAbundantIsotopeIndex] = new Tuple<Peak, int>(spec.Peaks[mostAbundantIsotopePeakIndex], mostAbundantIsotopePeakIndex);

            // go down
            var peakIndex = mostAbundantIsotopePeakIndex - 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i >= 0; i--)
                {
                    var peakMz = spec.Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        var peak = spec.Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Item1.Intensity)
                        {
                            observedPeaks[isotopeIndex] = new Tuple<Peak, int>(peak, peakIndex);
                        }
                    }
                }
            }

            // go up
            peakIndex = mostAbundantIsotopePeakIndex + 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i < spec.Peaks.Length; i++)
                {
                    var peakMz = spec.Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        var peak = spec.Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Item1.Intensity)
                        {
                            observedPeaks[isotopeIndex] = new Tuple<Peak, int>(peak, peakIndex);
                        }
                    }
                }
            }

            return observedPeaks;
        }

        // Select the best peak within +/- filteringWindowSize
        //public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
        //    Peak[] peaks, int minCharge, int maxCharge, 
        //    int isotopeOffsetTolerance, double filteringWindowSize,
        //    Tolerance tolerance, double corrScoreThreshold)
        //{

        //    var monoIsotopePeakList = new List<DeconvolutedPeak>();
        //    for (var peakIndex = 0; peakIndex < peaks.Length; peakIndex++)
        //    {
        //        var peak = peaks[peakIndex];

        //        // Check whether peak has the maximum intensity within the window
        //        var isBest = true;

        //        var prevIndex = peakIndex - 1;
        //        while (prevIndex >= 0)
        //        {
        //            var prevPeak = peaks[prevIndex];
        //            if ((peak.Mz - prevPeak.Mz) > filteringWindowSize) break;
        //            if (prevPeak.Intensity > peak.Intensity)
        //            {
        //                isBest = false;
        //                break;
        //            }
        //            prevIndex--;
        //        }

        //        if (!isBest) continue;

        //        var nextIndex = peakIndex + 1;
        //        while (nextIndex < peaks.Length)
        //        {
        //            var nextPeak = peaks[nextIndex];
        //            if ((nextPeak.Mz - peak.Mz) > filteringWindowSize) break;
        //            if (nextPeak.Intensity > peak.Intensity)
        //            {
        //                isBest = false;
        //                break;
        //            }
        //            nextIndex++;
        //        }

        //        if (!isBest) continue;

        //        // peak has the maximum intensity, window = [prevIndex+1,nextIndex-1]

        //        var window = new Peak[nextIndex - prevIndex - 1];
        //        Array.Copy(peaks, prevIndex + 1, window, 0, window.Length);
        //        var windowSpectrum = new Spectrum(window, 1);
        //        var peakMz = peak.Mz;

        //        double bestCorrelation = 0.0;
        //        DeconvolutedPeak bestPeak = null;

        //        for (var charge = maxCharge; charge >= minCharge; charge--)
        //        {
        //            var mass = (peak.Mz * charge) - charge * Constants.Proton;
        //            var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(mass);
        //            var mostAbundantIsotopeIndex = isotopomerEnvelope.MostAbundantIsotopeIndex;

        //            for (var isotopeIndex = mostAbundantIsotopeIndex - isotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex + isotopeOffsetTolerance; isotopeIndex++)
        //            {
        //                var monoIsotopeMass = Ion.GetMonoIsotopicMass(peakMz, charge, isotopeIndex);
        //                //var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(monoIsotopeMass);
        //                var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
        //                if (observedPeaks == null) continue;

        //                var envelop = isotopomerEnvelope.Envolope;
        //                var observedIntensities = new double[observedPeaks.Length];

        //                for (var i = 0; i < observedPeaks.Length; i++)
        //                {
        //                    var observedPeak = observedPeaks[i];
        //                    observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
        //                }

        //                var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelop, observedIntensities);
        //                var bcDist = sim.Item1;
        //                var corr = sim.Item2;

        //                //if (corr < corrScoreThreshold && bcDist > 0.03) continue;

        //                if (corr < corrScoreThreshold) continue;

        //                // monoIsotopeMass is valid
        //                if (corr >= bestCorrelation)
        //                {
        //                    bestCorrelation = corr;
        //                    bestPeak = new DeconvolutedPeak(monoIsotopeMass, observedIntensities[mostAbundantIsotopeIndex], charge, corr, bcDist, observedPeaks);
        //                }
        //            }
        //        }

        //        if (bestPeak != null)
        //        {
        //            monoIsotopePeakList.Add(bestPeak);
        //        }
        //    }

        //    monoIsotopePeakList.Sort();
        //    return monoIsotopePeakList;
        //}

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
