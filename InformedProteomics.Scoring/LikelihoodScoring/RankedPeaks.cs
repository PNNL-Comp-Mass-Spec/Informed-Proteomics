using System;
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
        public List<IonType> Iontypes { get; private set; }
        public RankInfo(double intensity)
        {
            Intensity = intensity;
            Iontypes = new List<IonType>();
        }
    }

    public class RankedPeaks
    {
        private readonly Spectrum _spectrum;
        public List<RankInfo> Ranks { get; private set; }
        public Dictionary<IonType, double> NotFound { get; private set; }
        public int NotFoundTotal { get; private set; }

        public RankedPeaks(Spectrum spectrum)
        {
            _spectrum = spectrum;
            var peaks = _spectrum.Peaks;
            Ranks = new List<RankInfo>();
            foreach (var peak in peaks)
            {
                Ranks.Add(new RankInfo(peak.Intensity));
            }
            Ranks.Add(new RankInfo(0));
            Ranks.Sort(new RankComparer());
            NotFound = new Dictionary<IonType, double>();
            NotFoundTotal = 0;
        }

        public int RankIon(IonType ionType, Ion ion, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance);
            int position;
            if (peak != null)
            {
                var intensity = peak.Intensity;
                var search = new RankInfo(intensity);
                position = Ranks.BinarySearch(search, new RankComparer());
                Ranks[position].Iontypes.Add(ionType);
                position++;
            }
            else
            {
                if (!NotFound.ContainsKey(ionType))
                    NotFound.Add(ionType, 0.0);
                NotFound[ionType]++;
                NotFoundTotal++;
                position = -1;
            }
            return position;
        }

        public int RankMz(double mz, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(mz, tolerance);
            if (peak == null) return -1;

            var intensity = peak.Intensity;
            var search = new RankInfo(intensity);
            var position = Ranks.BinarySearch(search, new RankComparer());
            return position;
        }
    }
}
