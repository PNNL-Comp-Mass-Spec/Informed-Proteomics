using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public interface ISpectrumExtractor
    {
        Spectrum GetSpectrum(int scanNum);
    }
}
