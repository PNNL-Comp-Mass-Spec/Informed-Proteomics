
namespace InformedProteomics.TopDown.Quantification
{
    public class MSDeconvNode
    {
        public MSDeconvNode(int scanNumber, double realMonoMass, double realIntensitySum, int charge)
        {
            ScanNumber = scanNumber;
            RealMonoMass = realMonoMass;
            RealIntensitySum = realIntensitySum;
            Charge = charge;
        }

        public int ScanNumber { get; }
        public int Charge { get; }
        public double RealMonoMass { get; }
        public double RealIntensitySum { get; }
    }
}
