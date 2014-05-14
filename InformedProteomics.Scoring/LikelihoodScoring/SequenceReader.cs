using System.Linq;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    class SequenceReader: ISequenceReader
    {
        public Sequence GetSequence(string sequence)
        {
            Sequence seq;
            if (sequence.Contains('('))
            {
                var sequenceReader = new MgfSequenceReader();
                seq = sequenceReader.GetSequence(sequence);
            }
            else
            {
                seq = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequence);
            }
            return seq;
        }
    }
}
