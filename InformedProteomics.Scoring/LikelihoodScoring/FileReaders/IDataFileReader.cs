using System.Collections.Generic;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    public interface IDataFileReader
    {
        IList<SpectrumMatch> Read();
    }
}
