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
        public Dictionary<IonType, double> NotFound { get; private set; }
        public int NotFoundTotal { get; private set; }

        public RankedPeaks(Spectrum spectrum)
        {
            _spectrum = spectrum;
            var peaks = _spectrum.Peaks;
            Ranks = new List<RankInfo>();
            foreach (var peak in peaks)
            {
                Ranks.Add(new RankInfo(peak.Intensity, null));
            }
            Ranks.Add(new RankInfo(0, null));
            Ranks.Sort(new RankComparer());
            NotFound = new Dictionary<IonType, double>();
            NotFoundTotal = 0;
        }

        public int RankIon(IonType ionType, Ion ion, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance);
            var position = 0;
            if (peak != null)
            {
                var intensity = peak.Intensity;
                var search = new RankInfo(intensity, null);
                position = Ranks.BinarySearch(search, new RankComparer());
            }
            else
            {
                if (!NotFound.ContainsKey(ionType))
                    NotFound.Add(ionType, 0.0);
                NotFound[ionType]++;
                NotFoundTotal++;
                return -1;
            }
            Ranks[position].Iontype = ionType;
            return position;
        }

        public int RankMz(double mz, Tolerance tolerance)
        {
            var peak = _spectrum.FindPeak(mz, tolerance);
            if (peak == null) return -1;

            var intensity = peak.Intensity;
            var search = new RankInfo(intensity, null);
            var position = Ranks.BinarySearch(search, new RankComparer());
            return position;
        }
    }
}
