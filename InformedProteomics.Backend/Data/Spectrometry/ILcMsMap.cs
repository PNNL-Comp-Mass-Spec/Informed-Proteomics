using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface ILcMsMap
    {
        // return list of (charge, MS1_ScanNumber)
        IEnumerable<List<Tuple<int, int>>> GetProbableChargeScanRegions(double monoIsotopicMass);
    }
}
