namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IDeconvolutedSpectrum
    {
        DeconvolutedPeak FindPeak(Composition.Composition composition, Tolerance tolerance);
    }
}
