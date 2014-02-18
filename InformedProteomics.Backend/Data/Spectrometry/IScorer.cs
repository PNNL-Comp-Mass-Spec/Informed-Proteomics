using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IScorer
    {
        double GetPrecursorIonScore(Ion precursorIon);
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition);
    }
}
