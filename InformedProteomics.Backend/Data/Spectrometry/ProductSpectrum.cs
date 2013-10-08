using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ProductSpectrum: Spectrum
    {
        public ProductSpectrum(double[] mzArr, double[] intensityArr, int scanNum) : base(mzArr, intensityArr, scanNum)
        {
        }

        public ProductSpectrum(ICollection<Peak> peaks, int scanNum): base(peaks, scanNum)
        {
        }

        public ActivationMethod ActivationMethod { get; set; }
        public IsolationInfo IsolationInfo { get; set; }

        public void SetMsLevel(int msLevel)
        {
            MsLevel = msLevel;
        }
    }
}
