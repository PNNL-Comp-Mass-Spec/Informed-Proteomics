using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Summd spectrum
    /// </summary>
    public class SummedSpectrum : Spectrum
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="scanNum"></param>
        public SummedSpectrum(ICollection<Peak> peaks, int scanNum) : base(peaks, scanNum)
        {
        }

        /// <summary>
        /// Scan numbers that were summed to make this <see cref="SummedSpectrum"/>
        /// </summary>
        public IList<int> ScanNums { get; set; }
    }
}
