using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// A fragment of a sequence
    /// </summary>
    [Obsolete("This class appears to be unused")]
    public class Fragment
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ionSymbol"></param>
        public Fragment(string ionSymbol)
        {
            IonSymbol = ionSymbol;
        }

        /// <summary>
        /// The charge state of the fragment
        /// </summary>
        public int ChargeState { get; set; }

        /// <summary>
        /// The m/z of the fragment
        /// </summary>
        public double Mz { get; set; }

        /// <summary>
        /// The mass of the fragment
        /// </summary>
        public double Mass { get; set; }

        /// <summary>
        /// If this fragment is y6, the FragmentIonClassBase is y
        /// </summary>
        public string IonType { get; set; }

        /// <summary>
        /// The index of the residue this fragment breaks the peptide at
        /// </summary>
        public int ResidueNumber { get; set; }

        /// <summary>
        /// The full ion symbol of the fragment; e.g. y6 or y6++ or y6-H2O++
        /// </summary>
        public string IonSymbol { get; }

        /// <summary>
        /// Check if two fragments are equal
        /// </summary>
        /// <param name="other"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(Fragment other)
        {
            if (other == null)
            {
                return false;
            }

            if (ReferenceEquals(this, other))
            {
                return true;
            }

            return Equals(other.IonSymbol, IonSymbol);
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (obj == null)
            {
                return false;
            }

            if (ReferenceEquals(this, obj))
            {
                return true;
            }

            if (obj.GetType() != typeof(Fragment))
            {
                return false;
            }

            return Equals((Fragment)obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return IonSymbol.GetHashCode();
        }

        /// <summary>
        /// Overloaded equality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        public static bool operator ==(Fragment left, Fragment right)
        {
            return Equals(left, right);
        }

        /// <summary>
        /// Overloaded inequality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        public static bool operator !=(Fragment left, Fragment right)
        {
            return !Equals(left, right);
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return string.Format("ChargeState: {0}, Mz: {1}, AveragineMass: {2}, ResidueNumber: {3}, FragmentIonClassBase: {4}, IonSymbol: {5}", ChargeState, Mz, Mass, ResidueNumber, IonType, IonSymbol);
        }
    }
}
