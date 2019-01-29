using System;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// MS2 spectrum isolation window data
    /// </summary>
    public class IsolationWindow : IComparable<IsolationWindow>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="isolationWindowTargetMz"></param>
        /// <param name="isolationWindowLowerOffset"></param>
        /// <param name="isolationWindowUpperOffset"></param>
        /// <param name="monoisotopicMz"></param>
        /// <param name="charge"></param>
        public IsolationWindow(
            double isolationWindowTargetMz,
            double isolationWindowLowerOffset,
            double isolationWindowUpperOffset,
            double? monoisotopicMz = null,
            int? charge = null
            )
        {
            IsolationWindowTargetMz = isolationWindowTargetMz;
            IsolationWindowLowerOffset = isolationWindowLowerOffset;
            IsolationWindowUpperOffset = isolationWindowUpperOffset;
            MonoisotopicMz = monoisotopicMz;
            Charge = charge;
        }

        /// <summary>
        /// Target m/z
        /// </summary>
        public double IsolationWindowTargetMz { get; }

        /// <summary>
        /// Lower mass offset
        /// </summary>
        public double IsolationWindowLowerOffset { get; }

        /// <summary>
        /// Upper mass offset
        /// </summary>
        public double IsolationWindowUpperOffset { get; }

        /// <summary>
        /// Precursor Monoisotopic m/z
        /// </summary>
        public double? MonoisotopicMz { get; set; }

        /// <summary>
        /// Precursor Charge
        /// </summary>
        public int? Charge { get; set; }

        /// <summary>
        /// Lower mass
        /// </summary>
        public double MinMz => IsolationWindowTargetMz - IsolationWindowLowerOffset;

        /// <summary>
        /// Upper mass
        /// </summary>
        public double MaxMz => IsolationWindowTargetMz + IsolationWindowUpperOffset;

        /// <summary>
        /// Isolation window width
        /// </summary>
        public double Width => IsolationWindowUpperOffset + IsolationWindowLowerOffset;

        /// <summary>
        /// Monoisotopic mass
        /// </summary>
        public double? MonoisotopicMass
        {
            get
            {
                if (MonoisotopicMz != null && Charge != null)
                    return (MonoisotopicMz - Constants.Proton) * Charge;

                return null;
            }
        }

        /// <summary>
        /// True of the mz is within the isolation window
        /// </summary>
        /// <param name="mz"></param>
        /// <returns></returns>
        public bool Contains(double mz)
        {
            return mz >= MinMz && mz < MaxMz;
        }

        /// <summary>
        /// Object equality
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        protected bool Equals(IsolationWindow other)
        {
            if (Math.Abs(MinMz - other.MinMz) < 0.01 && Math.Abs(MaxMz - other.MaxMz) < 0.01) return true;
            return false;
        }

        /// <summary>
        /// Comparer for sorting, orders by Target m/z
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(IsolationWindow other)
        {
            return IsolationWindowTargetMz.CompareTo(other.IsolationWindowTargetMz);
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (obj == null) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((IsolationWindow)obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            unchecked
            {
                var hashCode = IsolationWindowTargetMz.GetHashCode();
                hashCode = (hashCode * 397) ^ IsolationWindowLowerOffset.GetHashCode();
                hashCode = (hashCode * 397) ^ IsolationWindowUpperOffset.GetHashCode();
                return hashCode;
            }
        }
    }
}
