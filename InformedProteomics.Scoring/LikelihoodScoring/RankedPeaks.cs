using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    class RankComparer : IComparer<RankInfo>
    {
        int IComparer<RankInfo>.Compare(RankInfo l, RankInfo r)
        {
            return (r.Intensity.CompareTo(l.Intensity));
        }
    }

    public class RankInfo
    {
        public double Intensity { get; private set; }
        public IonType Iontype { get; set; }
        public RankInfo(double intensity, IonType ionType)
        {
            Intensity = intensity;
            Iontype = ionType;
        }
    }

    public class RankedPeaks
    {
        private readonly Spectrum _spectrum;
        public List<RankInfo> Ranks { get; private set; }

        public RankedPeaks(Spectrum spectrum)
        {
            _spectrum = spectrum;
            var peaks = _spectrum.Peaks;
            Ranks = new List<RankInfo>();
            foreach (var peak in peaks)
            {
                Ranks.Add(new RankInfo(peak.Intensity, null));
            }
            Ranks.Sort(new RankComparer());
        }

        public void RankIon(IonType ionType, Ion ion, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance);
            if (peak != null)
            {
                var intensity = peak.Intensity;
                var search = new RankInfo(intensity, null);
                var position = Ranks.BinarySearch(search, new RankComparer());
                Ranks[position].Iontype = ionType;
            }
        }
    }
}
