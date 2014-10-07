using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface ILcMsMap
    {
        IEnumerable<ChargeScanRange> GetProbableChargeScanRegions(double monoIsotopicMass);
    }
}
