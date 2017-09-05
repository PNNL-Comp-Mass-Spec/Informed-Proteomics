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
            _enzyme = enzyme;
            _ntt = ntt;
        }

        /// <summary>
        /// Enzyme
        /// </summary>
        public Enzyme Enzyme
        {
            get { return _enzyme; }
        }

        /// <summary>
        /// Number of tolerable termini
        /// </summary>
        public int Ntt
        {
            get { return _ntt; }
        }

        private readonly Enzyme _enzyme;
        private readonly int _ntt;
    }
}
