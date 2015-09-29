using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class DeconvolutedSpectrum : Spectrum
    {
        public DeconvolutedSpectrum(Spectrum originalSpec, DeconvolutedPeak[] peaks)
            : base(originalSpec.ScanNum)
        {
            Peaks = new DeconvolutedPeak[peaks.Length];
            peaks.CopyTo(Peaks, 0);
            MsLevel = originalSpec.MsLevel;

            var ms2Spec = originalSpec as ProductSpectrum;

            if (ms2Spec == null)
            {
                ActivationMethod = ActivationMethod.Unknown;
            }
            else
            {
                ActivationMethod = ms2Spec.ActivationMethod;
            }
            MsLevel = originalSpec.MsLevel;
        }

        public ActivationMethod ActivationMethod { get; private set; }
    }
}
