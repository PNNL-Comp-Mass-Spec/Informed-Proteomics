using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.BottomUp
{
    using InformedProteomics.Backend.Data.Sequence;

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
