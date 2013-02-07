using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Scoring
{
    public interface IScorer
    {
        float Score { get; }
        Sequence Seq { get; }
        DatabaseMultipleSubTargetResult MatchedResult { get; }

    }
}
