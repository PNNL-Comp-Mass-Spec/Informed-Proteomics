using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumFilter
    {
        /// <summary>
        /// For each peak in the spectrum, get all peaks within mzWidth m/z of the peak.
        /// Keep the peak if it is one of the top retentionCount peaks in that range.
        /// </summary>
        /// <param name="spectrum">The spectrum of peaks to filter.</param>
        /// <param name="mzWidth">The the m/z distance to look at in each direction around each peak.</param>
        /// <param name="retentionCount">The number of peaks to keep within the range [peakMz-mzWidth, peakMz+mzWidth]</param>
        /// <returns>Spectrum of filtered peaks.</returns>
        public static Spectrum FilterNoisePeaks(Spectrum spectrum, int mzWidth=100, int retentionCount=6)
        {
            var peaks = spectrum.Peaks;
            var filteredPeaks = new List<Peak>();

            foreach (var peak in peaks)
            {
                var minMz = peak.Mz - mzWidth;
                var maxMz = peak.Mz + mzWidth;
                var searchRange = (from index in peaks
                                   where (index.Mz >= minMz && index.Mz <= maxMz)
                                   select index).ToList();
                if (searchRange.Count > retentionCount)
                {
                    searchRange.Sort(new ComparePeakByIntensity());
                    var topRange = searchRange.GetRange(0, retentionCount);
                    if (topRange.Contains(peak))
                        filteredPeaks.Add(peak);
                }
                else
                {
                    filteredPeaks.Add(peak);
                }
            }
            filteredPeaks.Sort();
            var filteredSpectrum = new Spectrum(filteredPeaks.ToArray(), spectrum.ScanNum);
            return filteredSpectrum;
        }
    }

    class ComparePeakByIntensity : IComparer<Peak>
    {
        public int Compare(Peak x, Peak y)
        {
            return x.Intensity.CompareTo(y.Intensity);
        }
    }
}
