using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using InformedProteomics.Backend.Data.Sequence;

    public interface IScorer
    {
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition, AminoAcid nTerminalResidue = null, AminoAcid cTerminalResidue = null);
    }
}
