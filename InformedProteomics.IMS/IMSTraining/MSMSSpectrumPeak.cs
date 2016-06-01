using System;
using System.Collections.Generic;

namespace InformedProteomics.IMS.IMSTraining
{
    public class MSMSSpectrumPeak : IComparable<MSMSSpectrumPeak>
    {
        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public MSMSSpectrumPeak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public int CompareTo(MSMSSpectrumPeak other)
        {
            if (Mz > other.Mz) return 1;
            if (other.Mz > Mz) return -1;

            if (Intensity > other.Intensity) return 1;
            if (other.Intensity > Intensity) return -1;
            return 0;
        }

        internal class IntensityComparer : IComparer<MSMSSpectrumPeak>
        {
            public int Compare(MSMSSpectrumPeak x, MSMSSpectrumPeak y)
            {
                if (x.Intensity > y.Intensity) return 1;
                if (y.Intensity > x.Intensity) return -1;

                if (x.Mz > y.Mz) return 1;
                if (x.Mz < y.Mz) return -1;
                return 0;
            }
        }
    }
}
