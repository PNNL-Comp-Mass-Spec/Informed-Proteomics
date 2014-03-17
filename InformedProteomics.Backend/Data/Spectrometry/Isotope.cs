using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Isotope: IComparable<Isotope>
    {
        public Isotope(int index, double ratio)
        {
            Index = index;
            Ratio = ratio;
        }

        public int Index { get; private set; }
        public double Ratio { get; private set; }
        public int CompareTo(Isotope other)
        {
            return Index.CompareTo(other.Index);
        }
    }
}
