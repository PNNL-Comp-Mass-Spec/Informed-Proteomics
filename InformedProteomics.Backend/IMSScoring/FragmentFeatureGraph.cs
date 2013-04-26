using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureGraph : Dictionary<FeatureNode, List<FeatureEdge>>
    {
        public double Score { get; private set; }
        public FragmentFeatureGraph(ImsDataCached imsData, PrecursorFeatureNode precursorNode, Feature precursorFeature, Ion precursorIon, Composition cutComposition, GroupParameter parameter)
        {
            Add(precursorNode, new List<FeatureEdge>());
            var fragmentNodes = GetFragmentNodes(imsData, precursorFeature, cutComposition, precursorIon, parameter);
            if (fragmentNodes.Count == 0) return;
            UpdateEdges(this, fragmentNodes); // from precursor to any of fragment nodes
            var primeNode = (FragmentFeatureNode)this[precursorNode].ElementAt(0).RNode;
            var targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, false);
            UpdateEdges(this, targetNodes); // to the nodes of different terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, false, false);
            UpdateEdges(this, targetNodes); // to the nodes of the same terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, true);
            UpdateEdges(this, targetNodes); // to the nodes of differently charged ions

            Score = GetScore(this);
        }
        
        static private double GetScore(FragmentFeatureGraph graph)
        {
            return graph.Keys.SelectMany(node => graph[node]).Sum(edge => edge.Score);
        }

        static private List<FragmentFeatureNode> GetTargetNodes(FragmentFeatureNode primeNode, IEnumerable<FragmentFeatureNode> nodes, bool diffTerminal, bool diffCharge) // if diffCharge is true, diffTerminal is ignored 
        {
            var targetNodes = new List<FragmentFeatureNode>();
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
            var maxWeight = double.NegativeInfinity;
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

        static private List<FragmentFeatureNode> GetFragmentNodes(ImsDataCached imsData, Feature precursorFeature, Composition cutComposition, Ion precursorIon, GroupParameter parameter)
        {
            var ionTypes = SubScoreFactory.GetIonTypes(parameter);

            var nodes = new List<FragmentFeatureNode>();
            var suffixComposition = precursorIon.Composition - cutComposition;

            foreach (var ionType in ionTypes)
            {
                var composition = ionType.IsPrefixIon ? cutComposition : suffixComposition;
                nodes.Add(new FragmentFeatureNode(IsotopomerFeatures.GetFramentIsotopomerFeatures(imsData, composition, ionType, precursorFeature), ionType, parameter));
            }
            return nodes;
        } 

    }
}
