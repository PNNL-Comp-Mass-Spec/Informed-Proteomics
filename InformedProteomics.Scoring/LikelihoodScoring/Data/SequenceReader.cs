using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    internal class SequenceReader: ISequenceReader
    {
        private readonly DataFileFormat _format;
        public SequenceReader(DataFileFormat format)
        {
            _format = format;
        }
        public Sequence GetSequence(string sequence)
        {
            Sequence seq;
            if (_format == DataFileFormat.Mgf)
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
