using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class DeconvolutedSpectrum : ProductSpectrum
    {
        public DeconvolutedSpectrum(Spectrum originalSpec, DeconvolutedPeak[] peaks)
            : base(originalSpec.ScanNum)
        {
            MsLevel = originalSpec.MsLevel;
            NativeId = originalSpec.NativeId;
            ElutionTime = originalSpec.ElutionTime;
            TotalIonCurrent = originalSpec.TotalIonCurrent;

            var ms2Spec = originalSpec as ProductSpectrum;

            if (ms2Spec == null)
            {
                ActivationMethod = ActivationMethod.Unknown;
            }
            else
            {
                ActivationMethod = ms2Spec.ActivationMethod;
                IsolationWindow = new IsolationWindow(ms2Spec.IsolationWindow.IsolationWindowTargetMz, ms2Spec.IsolationWindow.IsolationWindowLowerOffset, ms2Spec.IsolationWindow.IsolationWindowUpperOffset, ms2Spec.IsolationWindow.MonoisotopicMz, ms2Spec.IsolationWindow.Charge);
            }

            var dPeaks = new DeconvolutedPeak[peaks.Length];
            peaks.CopyTo(dPeaks, 0);
            Peaks = dPeaks;
        }

        public DeconvolutedSpectrum(ICollection<DeconvolutedPeak> peaks, int scanNum) : base(scanNum)
        {
            var dPeaks = new DeconvolutedPeak[peaks.Count];
            peaks.CopyTo(dPeaks, 0);
            Peaks = dPeaks;
        }
    }
}
