using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

        public int ScanNumber { get; private set; }
        public int Charge { get; private set; }
        public double RealMonoMass { get; private set; }
        public double RealIntensitySum { get; private set; }
    }
}
