using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Scoring
{
    public interface IFragmentScorer
    {
        double GetFragmentScorer(Ion peptideIon, Composition suffixFragmentComposition);
    }
}
