using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for objects that can supply spectra
    /// </summary>
    public interface ISpectrumExtractor
    {
        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        Spectrum GetSpectrum(int scanNum, bool includePeaks = true);
    }
}
