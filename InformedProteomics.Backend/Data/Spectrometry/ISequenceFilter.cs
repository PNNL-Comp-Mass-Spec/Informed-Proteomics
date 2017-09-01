using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Interface for Sequence filters
    /// </summary>
    public interface ISequenceFilter
    {
        /// <summary>
        /// Get MS2 scans that have a precursor mass that matches <paramref name="sequenceMass"/>
        /// </summary>
        /// <param name="sequenceMass"></param>
        /// <returns></returns>
        IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass);
    }
}
