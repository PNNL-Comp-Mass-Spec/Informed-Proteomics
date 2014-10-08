using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class SummedSpectrum: Spectrum
    {
        public SummedSpectrum(ICollection<Peak> peaks, int scanNum) : base(peaks, scanNum)
        {
        }

        public IList<int> ScanNums { get; set; }
    }
}
