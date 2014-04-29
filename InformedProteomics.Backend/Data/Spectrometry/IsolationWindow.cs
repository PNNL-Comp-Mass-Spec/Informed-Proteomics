using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class IsolationWindow
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
    }
}
