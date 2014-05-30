using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public interface ISequenceReader
    {
        Sequence GetSequence(string sequence);
    }
}
