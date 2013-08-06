namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class XicPeak
    {
        public XicPeak(int scanNum, double intensity)
        {
            ScanNum = scanNum;
            Intensity = intensity;
        }

        public int ScanNum { get; private set; }
        public double Intensity { get; private set; }
    }
}
