namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Peptide spectrum match data
    /// </summary>
    public class PeptideSpectrumMatch
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="spectrum"></param>
        public PeptideSpectrumMatch(Sequence.Sequence peptide, ProductSpectrum spectrum)
        {
            Peptide = peptide;
            Spectrum = spectrum;
        }

        /// <summary>
        /// Peptide sequence
        /// </summary>
        public Sequence.Sequence Peptide { get; }

        /// <summary>
        /// Spectrum that matches <see cref="Peptide"/>
        /// </summary>
        public ProductSpectrum Spectrum { get; }
    }
}
