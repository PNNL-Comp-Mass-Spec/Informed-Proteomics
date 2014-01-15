using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.GeneratingFunction;

namespace InformedProteomics.Scoring.BottomUp
{
    public class PeptideScoringGraph: IScoringGraph
    {
        public int GetNodeScore(int nodeIndex)
        {
            var cleavageMass = nodeIndex/Constants.RescalingConstant;
            throw new System.NotImplementedException();
        }

        public IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex)
        {
            throw new System.NotImplementedException();
        }

        public int GetNumNodes()
        {
            throw new System.NotImplementedException();
        }

        private readonly AminoAcidSet _aminoAcidSet;
        private readonly Enzyme _enzyme;
        private readonly RankScorer _scorer;

    }
}
