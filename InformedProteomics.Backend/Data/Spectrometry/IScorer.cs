using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IScorer
    {
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition, AminoAcid nTerminalResidue = null, AminoAcid cTerminalResidue = null);
    }
}
