using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Deconvoluter
    {
        // Select the best peak within +/- filteringWindowSize
        public static List<DeconvolutedPeak> GetDeconvolutedPeaks(
            Spectrum spec, int minCharge, int maxCharge, 
            int isotopeOffsetTolerance, double filteringWindowSize,
            Tolerance tolerance, double corrScoreThreshold, double bcDistThreshold = 1.0)
        {
            var peaks = spec.Peaks;

            var monoIsotopePeakList = new List<DeconvolutedPeak>();
            //            var massHash = new HashSet<int>();
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
                var windowSpectrum = new Spectrum(window, spec.ScanNum);
                var peakMz = peak.Mz;

                //var hasToSkip = new bool[maxCharge - minCharge + 1, isotopeOffsetTolerance * 2 + 1];
                for (var charge = maxCharge; charge >= minCharge; charge--)
                {
                    var mass = peak.Mz * charge;
                    var mostAbundantIsotopeIndex = Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex;

                    for (var isotopeIndex = mostAbundantIsotopeIndex - isotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex + isotopeOffsetTolerance; isotopeIndex++)
                    {
                        //if (
                        //    hasToSkip[
                        //        charge - minCharge, isotopeIndex - mostAbundantIsotopeIndex + isotopeOffsetTolerance])
                        //    continue;

                        //var monoIsotopeMass = (peakMz - Constants.Proton) * charge - isotopeIndex * Constants.C13MinusC12;
                        var monoIsotopeMass = Ion.GetMonoIsotopicMass(peakMz, charge, isotopeIndex);
                        var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(monoIsotopeMass);
                        var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, charge, isotopomerEnvelope,
                            tolerance, 0.1);
                        if (observedPeaks == null) continue;

                        var envelop = isotopomerEnvelope.Envolope;
                        var observedIntensities = new double[observedPeaks.Length];

                        for (var i = 0; i < observedPeaks.Length; i++)
                        {
                            var observedPeak = observedPeaks[i];
                            observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                        }
                        var corr = FitScoreCalculator.GetPearsonCorrelation(envelop, observedIntensities);
                        if (corr < corrScoreThreshold) continue;
                        
                        var bcDist = FitScoreCalculator.GetBhattacharyyaDistance(envelop, observedIntensities);
                        if (bcDist > bcDistThreshold) continue;
                        

                        // monoIsotopeMass is valid
                        monoIsotopePeakList.Add(new DeconvolutedPeak(peak, monoIsotopeMass, charge));
                        //foreach (var factorCharge in SimpleMath.GetFactors(charge))
                        //{
                        //    if (factorCharge < minCharge) continue;
                        //    hasToSkip[factorCharge - minCharge, isotopeIndex - mostAbundantIsotopeIndex + isotopeOffsetTolerance] = true;
                        //}
                    }
                }
            }

            monoIsotopePeakList.Sort();
            return monoIsotopePeakList;
        }

        public static List<DeconvolutedPeak> FilterOut(List<DeconvolutedPeak> peaks, double windowSize, int topRankCutoff)
        {
                // select each peak if it is top n within window (-window,+window) around it
                //Spectrum retSpec = (Spectrum)s.clone();
                //retSpec.clear();    // remove all peaks 
                var retSpec = new List<DeconvolutedPeak>();
                for (int peakIndex = 0; peakIndex < peaks.Count; peakIndex++)
                {
                    int rank = 1;

                    var thisPeak = peaks[peakIndex];
                    var thisMass = thisPeak.Mass;
                    var thisInten = thisPeak.Intensity;

                    // move left
                    int prevIndex = peakIndex - 1;
                    while (prevIndex >= 0)
                    {
                        var prevPeak = peaks[prevIndex];
                        if (thisMass - prevPeak.Mass > windowSize) break;
                        if (prevPeak.Intensity > thisInten) rank++;
                        prevIndex--;
                    }

                    // move right
                    int nextIndex = peakIndex + 1;
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
