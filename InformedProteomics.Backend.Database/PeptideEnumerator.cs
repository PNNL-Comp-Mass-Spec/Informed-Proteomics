using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Database
{
    public class PeptideEnumerator
    {
        public PeptideEnumerator(Enzyme enzyme, int ntt)
        {
            _enzyme = enzyme;
            _ntt = ntt;
        }

        public Enzyme Enzyme
        {
            get { return _enzyme; }
        }

        public int Ntt
        {
            get { return _ntt; }
        }

        private readonly Enzyme _enzyme;
        private readonly int _ntt;
    }
}
