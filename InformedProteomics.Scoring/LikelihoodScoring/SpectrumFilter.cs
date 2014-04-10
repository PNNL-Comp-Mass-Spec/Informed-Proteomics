using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumFilter
    {
        public static Spectrum GetFilteredSpectrum(Spectrum spectrum, double mzWidth=100.0, int retentionCount=6)
        {
            var peaks = spectrum.Peaks;
            var filteredPeaks = new List<Peak>();

            for (var peakIndex = 0; peakIndex < peaks.Count(); peakIndex++)
            {
                var rank = 1;
                var peak = peaks[peakIndex];

                var prevIndex = peakIndex - 1;
                while (prevIndex >= 0)
                {
                    var prevPeak = peaks[prevIndex];
                    if ((peak.Mz - prevPeak.Mz) > mzWidth) break;
                    if (prevPeak.Intensity > peak.Intensity) rank++;
                    prevIndex--;
                }

                var nextIndex = peakIndex + 1;
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
