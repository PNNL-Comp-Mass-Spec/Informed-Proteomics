using System;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class IsolationWindow: IComparable<IsolationWindow>
    {
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

        public double IsolationWindowTargetMz { get; private set; }
        public double IsolationWindowLowerOffset { get; private set; }
        public double IsolationWindowUpperOffset { get; private set; }
        public double? MonoisotopicMz { get; set; }
        public int? Charge { get; set; }

        public double MinMz
        {
            get
            {
                return IsolationWindowTargetMz - IsolationWindowLowerOffset;
            }
        }

        public double MaxMz
        {
            get
            {
                return IsolationWindowTargetMz + IsolationWindowUpperOffset;
            }
        }

        public double Width
        {
            get { return IsolationWindowUpperOffset + IsolationWindowLowerOffset; }
        }

        public double? MonoisotopicMass
        {
            get
            {
                if (MonoisotopicMz != null && Charge != null) return (MonoisotopicMz - Constants.Proton)*Charge;
                return null;
            }
        }

        public bool Contains(double mz)
        {
            return mz >= MinMz && mz < MaxMz;
        }

        protected bool Equals(IsolationWindow other)
        {
            if (Math.Abs(MinMz - other.MinMz) < 0.01 && Math.Abs(MaxMz - other.MaxMz) < 0.01) return true;
            return false;
        }

        public int CompareTo(IsolationWindow other)
        {
            return IsolationWindowTargetMz.CompareTo(other.IsolationWindowTargetMz);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((IsolationWindow)obj);
        }

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
