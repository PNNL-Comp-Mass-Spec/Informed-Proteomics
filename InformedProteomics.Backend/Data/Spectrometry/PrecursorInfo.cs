using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class PrecursorInfo
    {
        public PrecursorInfo(
            int precursorScan,
            double isolationWindowTargetMz,
            double isolationWindowLowerOffset,
            double isolationWindowUpperOffset
            )
        {
            PrecursorScan = precursorScan;
            IsolationWindowTargetMz = isolationWindowTargetMz;
            IsolationWindowLowerOffset = isolationWindowLowerOffset;
            IsolationWindowUpperOffset = isolationWindowUpperOffset;
        }

        public int PrecursorScan { get; private set; }
        public double IsolationWindowTargetMz { get; private set; }
        public double IsolationWindowLowerOffset { get; private set; }
        public double IsolationWindowUpperOffset { get; private set; }
    }
}
