using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MathAndStats;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Deconvolution class
    /// </summary>
    public class Deconvoluter
    {
        /// <summary>
        /// Largest mass that will be deconvoluted
        /// </summary>
        public const double MaxMass = 50000;

        private readonly int minCharge;

        private readonly int maxCharge;

        private readonly int isotopeOffsetTolerance;

        private readonly Tolerance tolerance;

        private readonly double corrScoreThreshold;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="isotopeOffsetTolerance"></param>
        /// <param name="tolerance"></param>
        /// <param name="corrScoreThreshold"></param>
        public Deconvoluter(
            int minCharge,
            int maxCharge,
            int isotopeOffsetTolerance,
            Tolerance tolerance,
            double corrScoreThreshold = 0.7)
        {
            this.minCharge = minCharge;
            this.maxCharge = maxCharge;
            this.isotopeOffsetTolerance = isotopeOffsetTolerance;
            this.tolerance = tolerance;
            this.corrScoreThreshold = corrScoreThreshold;
        }

        /// <summary>
        /// Get a deconvoluted spectrum that combines multiple charge states
        /// </summary>
        /// <param name="spectrum"></param>
        /// <returns></returns>
        public DeconvolutedSpectrum GetCombinedDeconvolutedSpectrum(Spectrum spectrum)
        {
            return GetCombinedDeconvolutedSpectrum(
                        spectrum,
                        minCharge,
                        maxCharge,
                        isotopeOffsetTolerance,
                        tolerance,
                        corrScoreThreshold);
        }

        /// <summary>
        /// Get the deconvoluted spectrum
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="isotopeOffsetTolerance"></param>
        /// <param name="filteringWindowSize"></param>
        /// <param name="tolerance"></param>
        /// <param name="corrScoreThreshold"></param>
        /// <returns></returns>
        public static DeconvolutedSpectrum GetDeconvolutedSpectrum(
                    Spectrum spec, int minCharge, int maxCharge,
                    int isotopeOffsetTolerance, double filteringWindowSize,
                    Tolerance tolerance, double corrScoreThreshold = 0.7)
        {
            var peaks = GetDeconvolutedPeaks(spec.ScanNum, spec.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize, tolerance, corrScoreThreshold);

            return new DeconvolutedSpectrum(spec, peaks.ToArray());
        }

        /// <summary>
        /// Get a deconvoluted spectrum that combines multiple charge states
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="isotopeOffsetTolerance"></param>
        /// <param name="tolerance"></param>
        /// <param name="corrScoreThreshold"></param>
        /// <returns></returns>
        public static DeconvolutedSpectrum GetCombinedDeconvolutedSpectrum(
            Spectrum spectrum,
            int minCharge,
            int maxCharge,
            int isotopeOffsetTolerance,
            Tolerance tolerance,
            double corrScoreThreshold)
        {
            spectrum.FilterNoiseByIntensityHistogram();
            var deconvolutedPeaks = GetDeconvolutedPeaks_new(spectrum.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, tolerance, corrScoreThreshold);
            var binHash = new Dictionary<int, DeconvolutedPeak>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = (int)Math.Round(mass * (2 + Constants.RescalingConstant));
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

        /// <summary>
        /// Get the deconvoluted peaks that correspond to the provided peak list
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="isotopeOffsetTolerance"></param>
        /// <param name="tolerance"></param>
        /// <param name="corrScoreThreshold"></param>
        /// <returns></returns>
        public static List<DeconvolutedPeak> GetDeconvolutedPeaks_new(
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

                var bestScore = 0.0;
                DeconvolutedPeak bestPeak = null;
                Tuple<Peak, int>[] bestObservedPeaks = null;

                for (var charge = minCharge; charge <= maxCharge; charge++)
                {
                    var mass = peak.Mz * charge - (charge * Constants.Proton);
                    if (mass > MaxMass)
                    {
                        continue;
                    }

                    var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(mass);
                    var mostAbundantIsotopeIndex = isotopomerEnvelope.MostAbundantIsotopeIndex;
                    var offsetTolerance = isotopeOffsetTolerance;
                    if (isotopeOffsetTolerance < 0)
                    {
                        offsetTolerance = isotopomerEnvelope.Envelope.Length;
                    }

                    for (var isotopeIndex = mostAbundantIsotopeIndex - offsetTolerance;
                         isotopeIndex <= mostAbundantIsotopeIndex + offsetTolerance;
                         isotopeIndex++)
                    {
                        var monoIsotopeMass = Ion.GetMonoIsotopicMass(peak.Mz, charge, isotopeIndex);

                        var observedPeaks = GetAllIsotopePeaks(spectrum, monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
                        if (observedPeaks == null)
                        {
                            continue;
                        }

                        var envelope = isotopomerEnvelope.Envelope;
                        var observedIntensities = new double[observedPeaks.Length];

                        var observedPeakCount = 0;
                        for (var i = 0; i < observedPeaks.Length; i++)
                        {
                            var observedPeak = observedPeaks[i];
                            if (observedPeak != null && peakUsed[observedPeak.Item2])
                            {
                                observedPeak = null;
                                observedPeaks[i] = null;
                            }

                            observedPeakCount += observedPeak != null ? 1 : 0;
                            observedIntensities[i] = observedPeak != null ? (float)observedPeak.Item1.Intensity : 0.0;
                        }

                        var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelope, observedIntensities);
                        var bcDist = sim.Item1;
                        var corr = sim.Item2;
                        var foundPeakRatio = observedPeakCount / ((double)envelope.Length);

                        var interferenceScore = 10.0;

                        var filteredObserved = observedPeaks.Where(p => p != null).ToArray();
                        if (filteredObserved.Length >= 2)
                        {
                            var allPeaks =
                                spectrum.Peaks.Where(p => p.Mz >= filteredObserved[0].Item1.Mz && p.Mz <= filteredObserved[filteredObserved.Length - 1].Item1.Mz).ToArray();
                            interferenceScore = CalculateInterferenceScore(allPeaks, filteredObserved);
                        }

                        bcDist = Math.Max(bcDist, double.Epsilon);

                        if (corr < corrScoreThreshold && bcDist > 0.1)
                        {
                            continue;
                        }

                        var score = (foundPeakRatio * corr) / (bcDist * (Math.Abs(mostAbundantIsotopeIndex - isotopeIndex) + 1) * interferenceScore);

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
                            bestPeak.ObservedPeakIndices.Add(p.Item2);
                            peakUsed[p.Item2] = true;
                        }
                    }
                }
            }

            monoIsotopePeakList.Sort();
            return monoIsotopePeakList;
        }

        private static double CalculateInterferenceScore(IReadOnlyList<Peak> peaks, IReadOnlyList<Tuple<Peak, int>> envelopePeaks, int numBins = 10)
        {
            if (envelopePeaks.Count < 2)
            {
                return 10.0;
            }

            var maxIntensity = envelopePeaks.Max(p => p.Item1.Intensity);

            var bins = new double[numBins];
            var peakIndex = 0;
            for (var envelopePeakIndex = 1; envelopePeakIndex < envelopePeaks.Count; envelopePeakIndex++)
            {
                var lowPeak = envelopePeaks[envelopePeakIndex - 1].Item1;
                var highPeak = envelopePeaks[envelopePeakIndex].Item1;
                var binSize = (highPeak.Mz - lowPeak.Mz) / numBins;
                var prevBin = 0;
                for (; peakIndex < peaks.Count; peakIndex++)
                {
                    var noisePeak = peaks[peakIndex];
                    if (noisePeak.Mz.Equals(lowPeak.Mz))
                    {
                        continue;
                    }

                    if (noisePeak.Mz.Equals(highPeak.Mz))
                    {
                        break;
                    }

                    for (var i = prevBin; i < numBins; i++)
                    {
                        var lowBinMz = lowPeak.Mz + binSize * i;
                        var highBinMz = lowPeak.Mz + binSize * (i + 1);
                        if ((i == 0 && noisePeak.Mz < lowBinMz) || (noisePeak.Mz >= lowBinMz && noisePeak.Mz <= highBinMz)
                            || (i == numBins - 1 && noisePeak.Mz > highBinMz))
                        {
                            prevBin = i;
                            bins[i] += noisePeak.Intensity / maxIntensity;
                        }
                    }
                }
            }

            var interferenceScore = 0.0;
            for (var i = 0; i < numBins; i++)
            {
                interferenceScore += Math.Pow(2, bins[i]);
            }

            return interferenceScore / envelopePeaks.Count;
        }

        private static Tuple<Peak, int>[] GetAllIsotopePeaks(
            Spectrum spec,
            double monoisotopicMass,
            int charge,
            IsotopomerEnvelope envelope,
            Tolerance tolerance,
            double relativeIntensityThreshold)
        {
            var mostAbundantIsotopeIndex = envelope.MostAbundantIsotopeIndex;
            var isotopomerEnvelope = envelope.Envelope;
            var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, mostAbundantIsotopeIndex);
            var mostAbundantIsotopePeakIndex = spec.FindPeakIndex(mostAbundantIsotopeMz, tolerance);
            if (mostAbundantIsotopePeakIndex < 0)
            {
                return null;
            }

            var observedPeaks = new Tuple<Peak, int>[isotopomerEnvelope.Length];
            observedPeaks[mostAbundantIsotopeIndex] = new Tuple<Peak, int>(spec.Peaks[mostAbundantIsotopePeakIndex], mostAbundantIsotopePeakIndex);

            // go down
            var peakIndex = mostAbundantIsotopePeakIndex - 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold)
                {
                    break;
                }

                var isotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsMz(isotopeMz);
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
                    if (peakMz <= maxMz)    // find match, move to previous isotope
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
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold)
                {
                    break;
                }

                var isotopeMz = Ion.GetIsotopeMz(monoisotopicMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsMz(isotopeMz);
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
                    if (peakMz >= minMz)    // find match, move to previous isotope
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

        /// <summary>
        /// Get the deconvoluted peaks, selecting the best peak within +/- filteringWindowSize
        /// </summary>
        /// <param name="scanNum">Scan number (included in any exceptions that are caught)</param>
        /// <param name="peaks"></param>
        /// <param name="minCharge"></param>
        /// <param name="maxCharge"></param>
        /// <param name="isotopeOffsetTolerance"></param>
        /// <param name="filteringWindowSize"></param>
        /// <param name="tolerance"></param>
        /// <param name="corrScoreThreshold"></param>
        /// <returns></returns>
        public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
            int scanNum, Peak[] peaks,
            int minCharge, int maxCharge,
            int isotopeOffsetTolerance, double filteringWindowSize,
            Tolerance tolerance, double corrScoreThreshold)
        {
            try
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
                        if ((peak.Mz - prevPeak.Mz) > filteringWindowSize)
                        {
                            break;
                        }

                        if (prevPeak.Intensity > peak.Intensity)
                        {
                            isBest = false;
                            break;
                        }
                        prevIndex--;
                    }

                    if (!isBest)
                    {
                        continue;
                    }

                    var nextIndex = peakIndex + 1;
                    while (nextIndex < peaks.Length)
                    {
                        var nextPeak = peaks[nextIndex];
                        if ((nextPeak.Mz - peak.Mz) > filteringWindowSize)
                        {
                            break;
                        }

                        if (nextPeak.Intensity > peak.Intensity)
                        {
                            isBest = false;
                            break;
                        }
                        nextIndex++;
                    }

                    if (!isBest)
                    {
                        continue;
                    }

                    // peak has the maximum intensity, window = [previousIndex+1, nextIndex-1]

                    var window = new Peak[nextIndex - prevIndex - 1];
                    Array.Copy(peaks, prevIndex + 1, window, 0, window.Length);
                    var windowSpectrum = new Spectrum(window, 1);
                    var peakMz = peak.Mz;

                    if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
                    {
                        var bestScore = 0.0;
                        DeconvolutedPeak bestPeak = null;

                        for (var charge = maxCharge; charge >= minCharge; charge--)
                        {
                            var mass = (peak.Mz * charge) - charge * Constants.Proton;
                            var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(mass);
                            var mostAbundantIsotopeIndex = isotopomerEnvelope.MostAbundantIsotopeIndex;

                            for (var isotopeIndex = mostAbundantIsotopeIndex - isotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex + isotopeOffsetTolerance; isotopeIndex++)
                            {
                                var monoIsotopeMass = Ion.GetMonoIsotopicMass(peakMz, charge, isotopeIndex);
                                var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
                                if (observedPeaks == null)
                                {
                                    continue;
                                }

                                var envelope = isotopomerEnvelope.Envelope;
                                var observedIntensities = new double[observedPeaks.Length];

                                for (var i = 0; i < observedPeaks.Length; i++)
                                {
                                    var observedPeak = observedPeaks[i];
                                    observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                                }

                                var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelope, observedIntensities);
                                var bcDist = sim.Item1;
                                var corr = sim.Item2;
                                var score = corr / (bcDist * ((double)Math.Abs(isotopeIndex - mostAbundantIsotopeIndex) / envelope.Length));

                                if (corr < corrScoreThreshold && bcDist > 0.03)
                                {
                                    continue;
                                }

                                // monoIsotopeMass is valid
                                if (score >= bestScore)
                                {
                                    bestScore = score;
                                    bestPeak = new DeconvolutedPeak(monoIsotopeMass, observedIntensities[mostAbundantIsotopeIndex], charge, corr, bcDist, observedPeaks);
                                }
                            }
                        }

                        if (bestPeak != null)
                        {
                            monoIsotopePeakList.Add(bestPeak);
                        }
                    }
                    else
                    {
                        for (var charge = maxCharge; charge >= minCharge; charge--)
                        {
                            var mass = (peak.Mz * charge) - charge * Constants.Proton;
                            var mostAbundantIsotopeIndex = Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex;

                            for (var isotopeIndex = mostAbundantIsotopeIndex - isotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex + isotopeOffsetTolerance; isotopeIndex++)
                            {
                                var monoIsotopeMass = Ion.GetMonoIsotopicMass(peakMz, charge, isotopeIndex);
                                var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(monoIsotopeMass);
                                var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, charge, isotopomerEnvelope, tolerance, 0.1);
                                if (observedPeaks == null)
                                {
                                    continue;
                                }

                                var envelope = isotopomerEnvelope.Envelope;
                                var observedIntensities = new double[observedPeaks.Length];

                                for (var i = 0; i < observedPeaks.Length; i++)
                                {
                                    var observedPeak = observedPeaks[i];
                                    observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                                }

                                var sim = FitScoreCalculator.GetDistanceAndCorrelation(envelope, observedIntensities);
                                var bcDist = sim.Item1;
                                var corr = sim.Item2;

                                if (corr < corrScoreThreshold && bcDist > 0.03)
                                {
                                    continue;
                                }

                                var deconvPeak = new DeconvolutedPeak(monoIsotopeMass, observedIntensities[mostAbundantIsotopeIndex], charge, corr, bcDist, observedPeaks);
                                monoIsotopePeakList.Add(deconvPeak);
                            }
                        }
                    }
                }

                monoIsotopePeakList.Sort();
                return monoIsotopePeakList;
            }
            catch (Exception ex)
            {
                throw new Exception(string.Format("Error getting deconvoluted peaks for scan {0} in GetDeconvolutedPeaks: {1}", scanNum, ex.Message), ex);
            }
        }

        /// <summary>
        /// Filters out some data based on parameters, maybe?
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="windowSize"></param>
        /// <param name="topRankCutoff"></param>
        /// <returns></returns>
        public static List<DeconvolutedPeak> FilterOut(List<DeconvolutedPeak> peaks, double windowSize, int topRankCutoff)
        {
            var retSpec = new List<DeconvolutedPeak>();
            for (var peakIndex = 0; peakIndex < peaks.Count; peakIndex++)
            {
                var rank = 1;
                var thisPeak = peaks[peakIndex];
                var thisMass = thisPeak.Mass;
                var thisIntensity = thisPeak.Intensity;

                // move left
                var prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if (thisMass - prevPeak.Mass > windowSize)
                    {
                        break;
                    }

                    if (prevPeak.Intensity > thisIntensity)
                    {
                        rank++;
                    }

                    prevIndex--;
                }

                // move right
                for (var nextIndex = peakIndex + 1; nextIndex < peaks.Count; nextIndex++)
                {
                    var nextPeak = peaks[nextIndex];
                    if (nextPeak.Mass - thisMass > windowSize)
                    {
                        break;
                    }

                    if (nextPeak.Intensity > thisIntensity)
                    {
                        rank++;
                    }
                }
                if (rank <= topRankCutoff)
                {
                    retSpec.Add(thisPeak);
                }
            }

            return retSpec;
        }
    }
}
