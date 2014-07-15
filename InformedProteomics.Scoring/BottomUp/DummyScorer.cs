using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.BottomUp
{
    public class DummyScorer: IScorer
    {
        public double GetPrecursorIonScore(Ion precursorIon)
        {
            return 0.0;
        }

        public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
        {
            return 1;
        }
    }
}
