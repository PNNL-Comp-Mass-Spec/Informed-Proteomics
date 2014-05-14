using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public interface ISequenceReader
    {
        Sequence GetSequence(string sequence);
    }
}
