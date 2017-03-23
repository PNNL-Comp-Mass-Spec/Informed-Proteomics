using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.BottomUp
{
    public class DummyScorer: IScorer
    {
        public double GetPrecursorIonScore(Ion precursorIon)
        {
            return 0.0;
        }

        public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
        {
            return 1;
        }
    }
}
