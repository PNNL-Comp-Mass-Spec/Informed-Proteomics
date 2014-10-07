using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ChargeScanPair
    {
        public int Charge { get; private set; }
        public int Ms1ScanNumber { get; private set; }

        internal ChargeScanPair(int charge, int ms1)
        {
            Charge = charge;
            Ms1ScanNumber = ms1;
        }
    }
    
    public interface ILcMsMap
    {
        IEnumerable<List<ChargeScanPair>> GetProbableChargeScanRegions(double monoIsotopicMass);
    }
}
