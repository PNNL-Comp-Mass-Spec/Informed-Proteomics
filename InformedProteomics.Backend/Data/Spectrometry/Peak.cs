using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Peak data
    /// </summary>
    public class Peak: IComparable<Peak>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="intensity"></param>
        public Peak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        /// <summary>
        /// Peak m/z
        /// </summary>
        public double Mz { get; }

        /// <summary>
        /// Peak intensity
        /// </summary>
        public double Intensity { get; }

        /// <summary>
        /// Compare 2 peaks (for sorting by m/z)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Peak other)
        {
            var t = this as LcMsPeak;
            var o = other as LcMsPeak;
            if (t != null && o != null)
            {
                return t.CompareTo(o);
            }

            var mzCompare = Mz.CompareTo(other.Mz);
            if (mzCompare == 0)
            {
                // Only to force a stable sort
                return Intensity.CompareTo(other.Intensity);
            }
            return mzCompare;
        }

        /// <summary>
        /// Test 2 peaks for equality
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Peak other)
        {
            if (other == null)
            {
                return false;
            }

            if (ReferenceEquals(this, other))
            {
                return true;
            }

            return Math.Abs(Mz - other.Mz) < 1e-9;
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

            if (obj.GetType() != typeof (Peak))
            {
                return false;
            }

            return Equals((Peak) obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return Mz.GetHashCode();
        }

        /// <summary>
        /// Overloaded equality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator ==(Peak left, Peak right)
        {
            return Equals(left, right);
        }

        /// <summary>
        /// Overloaded inequality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator !=(Peak left, Peak right)
        {
            return !Equals(left, right);
        }
    }
}
