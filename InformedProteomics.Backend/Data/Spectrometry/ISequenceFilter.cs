using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface ISequenceFilter
    {
        IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass);
    }
}
