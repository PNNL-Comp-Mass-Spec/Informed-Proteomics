using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class XicPoint: IComparable<XicPoint>
    {
        public XicPoint(int scanNum, double mz, double intensity)
        {
            ScanNum = scanNum;
            Mz = mz;
            Intensity = intensity;
        }

        public int ScanNum { get; private set; }
        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public int CompareTo(XicPoint other)
        {
            return ScanNum - other.ScanNum;
        }
    }
}
