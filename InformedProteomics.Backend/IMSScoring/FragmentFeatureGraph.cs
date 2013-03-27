using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureGraph : Dictionary<FeatureNode, List<FeatureEdge>>
    {
        public float Score { get; private set; }
        public FragmentFeatureGraph(ImsDataCached imsData, FeatureNode precursorNode, Sequence peptide, int cutNumber)
        {
            Add(precursorNode, new List<FeatureEdge>());
            var fragmentNodes = GetFragmentNodes(imsData, precursorNode, peptide, cutNumber);
            if (fragmentNodes.Count == 0) return;
            UpdateEdges(this, fragmentNodes); // from precursor to any of fragment nodes
            var primeNode = this[precursorNode].ElementAt(0).RNode;
            var targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, false);
            UpdateEdges(this, targetNodes); // to the nodes of different terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, false, false);
            UpdateEdges(this, targetNodes); // to the nodes of the same terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, true);
            UpdateEdges(this, targetNodes); // to the nodes of differently charged ions

            Score = GetScore(this);
        }
        
        static private float GetScore(FragmentFeatureGraph graph)
        {
            return graph.Keys.SelectMany(node => graph[node]).Sum(edge => edge.Score);
        }

        static private List<FeatureNode> GetTargetNodes(FeatureNode primeNode, IEnumerable<FeatureNode> nodes, bool diffTerminal, bool diffCharge) // if diffCharge is true, diffTerminal is ignored 
        {
            var targetNodes = new List<FeatureNode>();
            //var primeNterm = primeNode.FragmentIonClassBase is IonType.NtermIonType;
            var primeNterm = primeNode.FragmentIonClassBase.IsPrefixIon;    // Modified by Sangtae

            foreach (var node in nodes)
            {
                if (node.Equals(primeNode)) continue;
                if (diffCharge)
                {
                    if(node.FragmentIonClassBase.Charge != primeNode.FragmentIonClassBase.Charge)
                        targetNodes.Add(node);
                }else
                {
                    if (node.FragmentIonClassBase.Charge != primeNode.FragmentIonClassBase.Charge) continue;
                    if (diffTerminal && primeNterm || !diffTerminal && !primeNterm)
                    {
                        if(node.FragmentIonClassBase.IsPrefixIon)
                            targetNodes.Add(node);
                    }
                    else
                    {
                        if(node.FragmentIonClassBase.IsPrefixIon)
                            targetNodes.Add(node);
                    }
                }
            }
            return targetNodes;
        } 

        static private void UpdateEdges(FragmentFeatureGraph graph, IEnumerable<FeatureNode> nodes)
        {
            var maxWeight = float.NegativeInfinity;
            FeatureEdge edgeWithMaxWeight = null;
            foreach (var rnode in nodes)
            {
                foreach (var lnode in graph.Keys)
                {
                    var edge = new FeatureEdge(lnode, rnode);
                    if (maxWeight < edge.Weight)
                    {
                        maxWeight = edge.Weight;
                        edgeWithMaxWeight = edge;
                    }
                }
            }
            if (edgeWithMaxWeight != null) graph[edgeWithMaxWeight.LNode].Add(edgeWithMaxWeight); 
        }

        static private List<FeatureNode> GetFragmentNodes(ImsDataCached imsData, FeatureNode precursorNode, Sequence peptide, int cutNumber)
        {
            var parameter = new FragmentParameter(peptide, cutNumber);
            var ionTypes = SubScoreFactory.GetIonTypes(parameter, precursorNode.FragmentIonClassBase.Charge);

            var nodes = new List<FeatureNode>();
            var prefixComposition = peptide.GetComposition(0, cutNumber);
            var suffixComposition = peptide.Composition - prefixComposition;

            foreach (var ionType in ionTypes)
            {
                double mz;
                if (ionType.IsPrefixIon)
                    mz = ionType.GetMz(prefixComposition);
                else
                    mz = ionType.GetMz(suffixComposition);
                nodes.Add(new FeatureNode(imsData.GetFramentFeature(mz, precursorNode.Feature), ionType, parameter));
            }
            return nodes;
        } 

    }
}
