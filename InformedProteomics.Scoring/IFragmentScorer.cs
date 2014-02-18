using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring
{
    public interface IFragmentScorer
    {
        double GetFragmentScorer(Ion peptideIon, Composition suffixFragmentComposition);
    }
}
