using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface ILcMsMap
    {
        IEnumerable<Ms1Feature> GetProbableChargeScanRegions(double monoIsotopicMass);
    }
}
