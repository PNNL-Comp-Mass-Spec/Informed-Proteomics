using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public class Ms1Feature
    {
        public Ms1Feature(int scanNum, double mass, int[] observedCharges, double abundance, ICollection<DeconvolutedPeak> peaks = null)
        {
            Mass = mass;
            Abundance = abundance;
            Charges = observedCharges;
            DeconvolutedPeaks = peaks;
            ScanNum = scanNum;
        }

        public readonly int ScanNum;
        public readonly int[] Charges;
        public double Mass { get; private set; }
        public int MinCharge { get { return Charges.Min(); } }
        public int MaxCharge { get { return Charges.Max(); } }
        public double Abundance { get; private set; }
        public ICollection<DeconvolutedPeak> DeconvolutedPeaks { get; private set; }
    }
}
