using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Isotope: IComparable<Isotope>
    {
        public Isotope(int index, float ratio)
        {
            Index = index;
            Ratio = ratio;
        }

        public int Index { get; private set; }
        public float Ratio { get; private set; }
        public int CompareTo(Isotope other)
        {
            return Index.CompareTo(other.Index);
        }
    }
}
