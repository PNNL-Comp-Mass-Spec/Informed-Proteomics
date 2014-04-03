using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumFilter
    {
        public static Spectrum FilterNoisePeaks(Spectrum spectrum, int mzWidth=100, int retentionCount=6)
        {
            var peaks = spectrum.Peaks;
            var filteredPeaks = new List<Peak>();

            var startingPeak = 0;

            while (startingPeak < peaks.Length)
            {
                var currentMz = peaks[startingPeak].Mz;
                var maxMz = currentMz + mzWidth;
                var peakSnapshot = new List<Peak>();
                int i;
                for (i = startingPeak; (i < peaks.Length && peaks[i].Mz <= maxMz); i++)
                {
                    peakSnapshot.Add(peaks[i]);
                }
                peakSnapshot.Sort(new ComparePeakByIntensity());
                int rangeWidth = retentionCount;
                if (peakSnapshot.Count < retentionCount)
                    rangeWidth = peakSnapshot.Count;
                filteredPeaks.AddRange(peakSnapshot.GetRange(0, rangeWidth));
                startingPeak = i;
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
