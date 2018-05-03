using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Class for hold MSn (product) spectrum information
    /// </summary>
    public class ProductSpectrum: Spectrum
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mzArr"></param>
        /// <param name="intensityArr"></param>
        /// <param name="scanNum"></param>
        public ProductSpectrum(IList<double> mzArr, IList<double> intensityArr, int scanNum) : base(mzArr, intensityArr, scanNum)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="scanNum"></param>
        public ProductSpectrum(ICollection<Peak> peaks, int scanNum): base(peaks, scanNum)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNum"></param>
        public ProductSpectrum(int scanNum) : base(scanNum)
        {
        }

        /// <summary>
        /// Activation Method
        /// </summary>
        public ActivationMethod ActivationMethod { get; set; }

        /// <summary>
        /// Isolation Window
        /// </summary>
        public IsolationWindow IsolationWindow { get; set; }

        /// <summary>
        /// Set the MS Level
        /// </summary>
        /// <param name="msLevel"></param>
        public void SetMsLevel(int msLevel)
        {
            MsLevel = msLevel;
        }
    }
}
