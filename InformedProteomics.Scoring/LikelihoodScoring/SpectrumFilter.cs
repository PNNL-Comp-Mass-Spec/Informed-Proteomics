using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumFilter
    {
        public static Spectrum FilterNoisePeaks(Spectrum spectrum, int mzWidth=100, int retentionCount=6)
        {
            var peaks = spectrum.Peaks;
            var filteredPeaks = new List<Peak>();

            for (int peakIndex = 0; peakIndex < peaks.Count(); peakIndex++)
            {
                int rank = 1;
                var peak = peaks[peakIndex];

                int prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if ((peak.Mz - prevPeak.Mz) > mzWidth) break;
                    if (prevPeak.Intensity > peak.Intensity) rank++;
                    prevIndex--;
                }

                int nextIndex = peakIndex + 1;
                while (nextIndex < peaks.Length)
                {
                    var nextPeak = peaks[nextIndex];
                    if ((nextPeak.Mz - peak.Mz) > mzWidth) break;
                    if (nextPeak.Intensity > peak.Intensity) rank++;
                    nextIndex++;
                }
                if (rank <= retentionCount) filteredPeaks.Add(peak);
            }
            var filteredSpectrum = new Spectrum(filteredPeaks.ToArray(), spectrum.ScanNum);
            return filteredSpectrum;
        }
    }
}
