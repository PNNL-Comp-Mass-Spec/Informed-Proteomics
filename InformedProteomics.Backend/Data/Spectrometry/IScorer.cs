using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IScorer
    {
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition);
    }
}
