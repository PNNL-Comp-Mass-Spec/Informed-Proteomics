using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IScorer
    {
        double GetPrecursorIonScore(Ion precursorIon);
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition);
    }
}
