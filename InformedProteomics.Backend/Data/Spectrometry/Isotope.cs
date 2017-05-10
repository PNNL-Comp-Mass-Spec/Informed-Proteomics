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

        public int Index { get; }

        public double Ratio { get; }

        /// <summary>
        /// Compare two isotopes by Index
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Isotope other)
        {
            return Index.CompareTo(other.Index);
        }

        public override string ToString()
        {
            return string.Format("{0} , {1:F4}", Index, Ratio);
        }
    }
}
