using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class PeptideSpectrumMatch
    {
        public PeptideSpectrumMatch(Sequence.Sequence peptide, ProductSpectrum spectrum)
        {
            Peptide = peptide;
            Spectrum = spectrum;
        }

        public Sequence.Sequence Peptide { get; private set; }
        public ProductSpectrum Spectrum { get; private set; }
    }
}
