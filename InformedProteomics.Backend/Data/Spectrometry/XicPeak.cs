using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class XicPeak: IComparable<XicPeak>
    {
        public XicPeak(int scanNum, double intensity)
        {
            ScanNum = scanNum;
            Intensity = intensity;
        }

        public int ScanNum { get; private set; }
        public double Intensity { get; private set; }
        public int CompareTo(XicPeak other)
        {
            return ScanNum - other.ScanNum;
        }
    }
}
