namespace InformedProteomics.Backend.Data.Composition
{
    /// <summary>
    /// Class for a single Isotopomer envelope
    /// </summary>
    public class IsotopomerEnvelope
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="envelope"></param>
        /// <param name="mostAbundantIsotopeIndex"></param>
        public IsotopomerEnvelope(double[] envelope, int mostAbundantIsotopeIndex)
        {
            Envelope = envelope;
            MostAbundantIsotopeIndex = mostAbundantIsotopeIndex;
        }

        /// <summary>
        /// The Isotopomer envelope
        /// </summary>
        public double[] Envelope { get; }

        /// <summary>
        /// Index in <see cref="Envelope"/> of the most abundant isotope
        /// </summary>
        public int MostAbundantIsotopeIndex { get; }
    }
}
