using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Execution
{
    public class SelectedMsMsFilter : ISequenceFilter
    {
        private readonly List<int> ms2Scans;

        public SelectedMsMsFilter(IEnumerable<int> ms2scans)
        {
            this.ms2Scans = new List<int>(ms2scans);
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return ms2Scans;
        }
    }
}
