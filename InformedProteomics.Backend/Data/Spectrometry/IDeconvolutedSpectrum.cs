namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Interface for deconvoluted spectra
    /// </summary>
    public interface IDeconvolutedSpectrum
    {
        /// <summary>
        /// Get the peak that matches <paramref name="composition"/> and <paramref name="tolerance"/>
        /// </summary>
        /// <param name="composition"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        DeconvolutedPeak FindPeak(Composition.Composition composition, Tolerance tolerance);
    }
}
