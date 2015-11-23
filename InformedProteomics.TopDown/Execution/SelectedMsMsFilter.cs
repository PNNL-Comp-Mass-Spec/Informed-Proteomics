using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.TopDown.Execution
{
    using InformedProteomics.Backend.Data.Spectrometry;
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
