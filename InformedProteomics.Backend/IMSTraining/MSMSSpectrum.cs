using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSTraining
{
    public class MSMSSpectrum : List<MSMSSpectrumPeak>
    {
        public Sequence Annotation { get; private set; }
        public int Charge { get; private set; }
        public double PrecursorMz { get; private set; }

        public MSMSSpectrum(int charge, double precursorMz)
        {
            Charge = charge;
            PrecursorMz = precursorMz;
        }

        public MSMSSpectrum(int charge, double precursorMz, Sequence annotation) : this(charge, precursorMz)
        {
            Annotation = annotation;
        }

        public MSMSSpectrumPeak GetPeak(double mz, Tolerance tolerance)
        {
            var matchList = GetPeaks(mz, tolerance);
            if (matchList.Count == 0) return null;

            var maxPeak = matchList[0];
            var intensityComparer = new MSMSSpectrumPeak.IntensityComparer();
            for (var i = 1; i < matchList.Count; i++)
            {
                var p = matchList[i];
                if (intensityComparer.Compare(maxPeak, p) < 0)
                    maxPeak = p;
            }
            return maxPeak;
        }

        public List<MSMSSpectrumPeak> GetPeaks(double mz, Tolerance tolerance)
        {
            var minMz = mz - tolerance.GetToleranceAsTh(mz);
            var maxMz = mz + tolerance.GetToleranceAsTh(mz);
            var matchList = new List<MSMSSpectrumPeak>();
            var start = BinarySearch(new MSMSSpectrumPeak(minMz, 0));
            if (start < 0) start = -start - 1;
            for (var i = start; i < Count; i++)
            {
                var p = this[i];
                if (p.Mz > maxMz) break;
                matchList.Add(p);
            }
            return matchList;
        }


    }
}
