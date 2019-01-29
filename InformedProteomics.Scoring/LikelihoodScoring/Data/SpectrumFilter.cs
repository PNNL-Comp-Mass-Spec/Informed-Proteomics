using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class SpectrumFilter
    {
        public static Spectrum GetFilteredSpectrum(Spectrum spectrum, double windowWidth=100.0, int retentionCount=6)
        {
            windowWidth = windowWidth/2;
            var peaks = spectrum.Peaks;
            var filteredPeaks = new List<Peak>();

            for (var peakIndex = 0; peakIndex < peaks.Length; peakIndex++)
            {
                var rank = 1;
                var peak = peaks[peakIndex];

                var prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if ((peak.Mz - prevPeak.Mz) > windowWidth) break;
                    if (prevPeak.Intensity > peak.Intensity) rank++;
                    prevIndex--;
                }

                var nextIndex = peakIndex + 1;
                while (nextIndex < peaks.Length)
                {
                    var nextPeak = peaks[nextIndex];
                    if ((nextPeak.Mz - peak.Mz) > windowWidth) break;
                    if (nextPeak.Intensity > peak.Intensity) rank++;
                    nextIndex++;
                }
                if (rank <= retentionCount) filteredPeaks.Add(peak);
            }
            var filteredSpectrum = new Spectrum(filteredPeaks.ToArray(), spectrum.ScanNum);
            return filteredSpectrum;
        }

        /// <summary>
        /// Filter out all peaks in a spectrum that are not explained by certain ion types.
        /// </summary>
        /// <param name="sequence">Sequence to calculate ions from.</param>
        /// <param name="spectrum">Spectrum to filter.</param>
        /// <param name="ionTypes">Ion types to find peaks for.</param>
        /// <param name="tolerance"></param>
        /// <returns>Filtered Peptide Spectrum Match</returns>
        public static Spectrum FilterIonPeaks(Sequence sequence, Spectrum spectrum, IonType[] ionTypes, Tolerance tolerance)
        {
            var filteredPeaks = new List<Peak>();
            var specMatch = new SpectrumMatch(sequence, spectrum);
            foreach (var ionType in ionTypes)
            {
                var ions = specMatch.GetCleavageIons(ionType);
                foreach (var ion in ions)
                {
                    var peak = spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance);
                    if (peak != null) filteredPeaks.Add(peak);
                }
            }
            filteredPeaks.Sort();
            return new Spectrum(filteredPeaks, spectrum.ScanNum) {MsLevel = 2};
        }
    }
}
