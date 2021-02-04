using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Deconvoluted spectrum
    /// </summary>
    public class DeconvolutedSpectrum : ProductSpectrum, IDeconvolutedSpectrum
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="originalSpec"></param>
        /// <param name="peaks"></param>
        public DeconvolutedSpectrum(Spectrum originalSpec, IReadOnlyList<DeconvolutedPeak> peaks)
            : base(originalSpec.ScanNum)
        {
            MsLevel = originalSpec.MsLevel;
            NativeId = originalSpec.NativeId;
            ElutionTime = originalSpec.ElutionTime;
            TotalIonCurrent = originalSpec.TotalIonCurrent;

            if (originalSpec is ProductSpectrum ms2Spec)
            {
                ActivationMethod = ms2Spec.ActivationMethod;
                IsolationWindow = new IsolationWindow(ms2Spec.IsolationWindow.IsolationWindowTargetMz,
                                                      ms2Spec.IsolationWindow.IsolationWindowLowerOffset,
                                                      ms2Spec.IsolationWindow.IsolationWindowUpperOffset, ms2Spec.IsolationWindow.MonoisotopicMz,
                                                      ms2Spec.IsolationWindow.Charge);
            }
            else
            {
                ActivationMethod = ActivationMethod.Unknown;
            }

            StorePeaks(peaks);
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peaks"></param>
        /// <param name="scanNum"></param>
        public DeconvolutedSpectrum(IEnumerable<DeconvolutedPeak> peaks, int scanNum) : base(scanNum)
        {
            StorePeaks(peaks.ToList());
        }

        /// <summary>
        /// Get the peak that matches <paramref name="composition"/> and <paramref name="tolerance"/>
        /// </summary>
        /// <param name="composition"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public DeconvolutedPeak FindPeak(Composition.Composition composition, Tolerance tolerance)
        {
            return base.FindPeak(composition.Mass, tolerance) as DeconvolutedPeak;
        }

        /// <summary>
        /// Store the peaks in the Peaks, assuring that each peak is an instance of DeconvolutedPeak
        /// </summary>
        /// <param name="peaks"></param>
        private void StorePeaks(IReadOnlyList<DeconvolutedPeak> peaks)
        {
            var dPeaks = new DeconvolutedPeak[peaks.Count];
            for (var i = 0; i < peaks.Count; i++)
            {
                dPeaks[i] = peaks[i];
            }

            // Store the deconvoluted peaks in Peaks
            // ReSharper disable once CoVariantArrayConversion
            Peaks = new DeconvolutedPeak[dPeaks.Length];

            for (var i = 0; i < dPeaks.Length; i++)
            {
                Peaks[i] = dPeaks[i];
            }
        }
    }
}
