using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureGraph
    {
        public List<FeatureNode> PrecursorFeatureNodes; // precursor feature nodes of different charge states
        public Sequence Peptide { get; private set; }
        private readonly Dictionary<FeatureNode, FragmentFeatureGraph>[] _fragmentFeatureGraphs; 

        public FeatureGraph(ImsDataCached imsData, List<FeatureNode> precursorFeatureNodes, Sequence peptide)
        {
            PrecursorFeatureNodes = precursorFeatureNodes;
            Peptide = peptide;

            _fragmentFeatureGraphs = new Dictionary<FeatureNode, FragmentFeatureGraph>[Peptide.Count];

            for (var cutNumber = 0; cutNumber < Peptide.Count; cutNumber++)
            {
                _fragmentFeatureGraphs[cutNumber] = new Dictionary<FeatureNode, FragmentFeatureGraph>();
                foreach (var precursorFeatureNode in precursorFeatureNodes)
                {
                    _fragmentFeatureGraphs[cutNumber][precursorFeatureNode] = new FragmentFeatureGraph(imsData, precursorFeatureNode, peptide, cutNumber);
                }
            }
        }

        public FragmentFeatureGraph GetFragmentFeatureGraph(FeatureNode precursorFeatureNode, int cutNumber)
        {
            return _fragmentFeatureGraphs[cutNumber][precursorFeatureNode];
        } 

    }
}
