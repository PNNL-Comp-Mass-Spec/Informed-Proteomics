using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    class SequenceReader: ISequenceReader
    {
        private readonly string _format;
        public SequenceReader(string format)
        {
            _format = format;
        }
        public Sequence GetSequence(string sequence)
        {
            Sequence seq;
            if (_format == "mgf")
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
