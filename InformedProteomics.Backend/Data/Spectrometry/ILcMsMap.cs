using System.Collections.Generic;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface ILcMsMap
    {
        IEnumerable<ChargeLcScanCluster> GetProbableChargeScanRegions(double monoIsotopicMass);
    }
}
