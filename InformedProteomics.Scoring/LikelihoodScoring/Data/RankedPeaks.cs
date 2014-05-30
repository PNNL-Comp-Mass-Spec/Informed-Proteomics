using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class RankedPeaks : Spectrum
    {
        public RankedPeaks(Spectrum spectrum)
            : base(spectrum.Peaks, spectrum.ScanNum)
        {
            _spectrum = spectrum;
            Array.Sort(Peaks, new IntensityComparer());
        }

        /// <summary>
        /// Find intensity rank of ion in spectrum.
        /// </summary>
        /// <param name="ion">Ion to search for.</param>
        /// <param name="tolerance"></param>
        /// <returns>Intensity rank of ion. 1-based.</returns>
        public int RankIon(Ion ion, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance);
            int position;
            if (peak != null)
            {
                var searchPeak = new Peak(peak.Mz, peak.Intensity);
                position = Array.BinarySearch(Peaks, searchPeak, new IntensityComparer());
                position++;
            }
            else
            {
                position = -1;
            }
            return position;
        }

        private readonly Spectrum _spectrum;
    }
}
