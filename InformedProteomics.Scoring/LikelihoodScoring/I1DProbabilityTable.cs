using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public interface I1DProbabilityTable<T1>
    {
        Probability<T1>[] GetProbabilities();
        T1[] GetBinEdges();
        void AddMatches(List<SpectrumMatch> matches);
    }
}
