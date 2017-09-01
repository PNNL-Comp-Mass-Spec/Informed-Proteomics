using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Isotope
    /// </summary>
    public class Isotope: IComparable<Isotope>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="index"></param>
        /// <param name="ratio"></param>
        public Isotope(int index, double ratio)
        {
            Index = index;
            Ratio = ratio;
        }

        /// <summary>
        /// Isotope index
        /// </summary>
        public int Index { get; }

        /// <summary>
        /// Isotope ratio
        /// </summary>
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

        /// <inheritdoc />
        public override string ToString()
        {
            return string.Format("{0} , {1:F4}", Index, Ratio);
        }
    }
}
