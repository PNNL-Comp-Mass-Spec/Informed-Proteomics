using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class XicPeak: IComparable<XicPeak>
    {
        public XicPeak(int scanNum, double mz, double intensity)
        {
            ScanNum = scanNum;
            Mz = mz;
            Intensity = intensity;
        }

        public int ScanNum { get; private set; }
        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public int CompareTo(XicPeak other)
        {
            return ScanNum - other.ScanNum;
        }
    }
}
