namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class PeptideSpectrumMatch
    {
        public PeptideSpectrumMatch(Sequence.Sequence peptide, ProductSpectrum spectrum)
        {
            Peptide = peptide;
            Spectrum = spectrum;
        }

        public Sequence.Sequence Peptide { get; }
        public ProductSpectrum Spectrum { get; }
    }
}
