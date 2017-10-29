using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Enumerate over peptides for bottom-up searches
    /// </summary>
    public class PeptideEnumerator
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="enzyme"></param>
        /// <param name="ntt"></param>
        public PeptideEnumerator(Enzyme enzyme, int ntt)
        {
            Enzyme = enzyme;
            Ntt = ntt;
        }

        /// <summary>
        /// Enzyme
        /// </summary>
        public Enzyme Enzyme { get; }

        /// <summary>
        /// Number of tolerable termini
        /// </summary>
        public int Ntt { get; }
    }
}
