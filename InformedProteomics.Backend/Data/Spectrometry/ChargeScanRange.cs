namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeScanRange
    {
        public ChargeScanRange(int minCharge, int maxCharge, int minScanNum, int maxScanNum)
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScanNum;
            MaxScanNum = maxScanNum;
        }

        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        public int MinScanNum { get; private set; }
        public int MaxScanNum { get; private set; }
    }
}
