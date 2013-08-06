using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Peak: IComparable<Peak>
    {
        public Peak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public int CompareTo(Peak other)
        {
            if (Mz < other.Mz) return -1;
            return Mz > other.Mz ? 1 : 0;
        }
    }

    /// <summary>
    /// Sort by reverse order of intensities (highest intensity comes first)
    /// </summary>
    public class IntensityComparer : IComparer<Peak>
    {
        public int Compare(Peak x, Peak y)
        {
            if (x.Intensity < y.Intensity) return 1;
            return x.Intensity > y.Intensity ? -1 : 0;
        }
    }
}
