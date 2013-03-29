using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.IMSScoring
{
    class FeatureGraphScorer
    {
        public float PrecursorFeatureNodesScore { get; private set; }
        public float[,] FragmentScore { get; private set; } // precursor node, cut number
        public float Score { get; private set; }
       
        public FeatureGraphScorer(FeatureGraph graph)
        {
            Score = PrecursorFeatureNodesScore = GetPrecursorFeatureNodesScore(graph.PrecursorFeatureNodes);
            FragmentScore = new float[graph.PrecursorFeatureNodes.Count,graph.Peptide.Count];
            for (var  i=0;i<graph.PrecursorFeatureNodes.Count;i++)
            {
                var precursorFeatureNode = graph.PrecursorFeatureNodes.ElementAt(i);
                for (var cutNumber = 0; cutNumber < graph.Peptide.Count; cutNumber++)
                {
                    FragmentScore[i, cutNumber] = graph.GetFragmentFeatureGraph(precursorFeatureNode, cutNumber).Score;
                    Score += FragmentScore[i, cutNumber];
                }
            } 
        }

        static private float GetPrecursorFeatureNodesScore(List<PrecursorFeatureNode> nodes)
        {
            return 0;//TODO
        }
    }
}
