namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeScanRange
    {
        public ChargeScanRange(int minCharge, int maxCharge, int minScanNum, int maxScanNum, int repScanNum, double score)
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScanNum;
            MaxScanNum = maxScanNum;
            RepresentativeScanNum = repScanNum;
            Score = score;
        }

        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        public int MinScanNum { get; private set; }
        public int MaxScanNum { get; private set; }

        public int RepresentativeScanNum { get; private set; }
        public double Score { get; private set; }
    }
}
