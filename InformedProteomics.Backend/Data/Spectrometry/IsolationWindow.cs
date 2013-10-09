using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class IsolationWindow
    {
        public IsolationWindow(
            double isolationWindowTargetMz,
            double isolationWindowLowerOffset,
            double isolationWindowUpperOffset
            )
        {
            IsolationWindowTargetMz = isolationWindowTargetMz;
            IsolationWindowLowerOffset = isolationWindowLowerOffset;
            IsolationWindowUpperOffset = isolationWindowUpperOffset;
        }

        public double IsolationWindowTargetMz { get; private set; }
        public double IsolationWindowLowerOffset { get; private set; }
        public double IsolationWindowUpperOffset { get; private set; }

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

        public bool Contains(double mz)
        {
            return mz >= MinMz && mz < MaxMz;
        }
    }
}
