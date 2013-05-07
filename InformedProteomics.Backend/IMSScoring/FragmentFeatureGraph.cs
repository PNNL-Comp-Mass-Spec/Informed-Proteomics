using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureGraph : Dictionary<FeatureNode, List<FeatureEdge>>
    {
        public double Score { get; private set; }
        public double NodeScore { get; private set; }
        public double RatioScore { get; private set; }
        public double LcScore { get; private set; }
        public double ImsScore { get; private set; }

        private readonly SubScoreFactory _scoringParams;

        public FragmentFeatureGraph(ImsDataCached imsData, PrecursorFeatureNode precursorNode, Feature precursorFeature, 
            Ion precursorIon, Composition cutComposition, GroupParameter parameter, SubScoreFactory scoringParams)
        {
            Add(precursorNode, new List<FeatureEdge>());
            var fragmentNodes = GetFragmentNodes(imsData, precursorFeature, cutComposition, precursorIon, parameter);

            //var nn = fragmentNodes.Count(no => no.Feature != null);
            //Console.WriteLine(this + "Num features : " + nn);

            _scoringParams = scoringParams;

            var usedNodes = new List<FeatureNode> {precursorNode};
            UpdateEdges(fragmentNodes, usedNodes); // from precursor to any of fragment nodes
            if (this[precursorNode].Count == 0)
            {
                NodeScore = RatioScore = LcScore = ImsScore = -1;
                Score = NodeScore + RatioScore + LcScore + ImsScore;
                return;
            }
            var primeNode = (FragmentFeatureNode) this[precursorNode][0].RNode;
            var targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, false);
            UpdateEdges(targetNodes, usedNodes); // to the nodes of different terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, false, false);
            UpdateEdges(targetNodes, usedNodes); // to the nodes of the same terminal ions
            targetNodes = GetTargetNodes(primeNode, fragmentNodes, true, true);
            UpdateEdges(targetNodes, usedNodes); // to the nodes of differently charged ions

            GetScore();
        }
        
        private void GetScore()
        {
            NodeScore = RatioScore = LcScore = ImsScore = 0;
            foreach (var node in Keys)
            {
                foreach (var edge in this[node])
                {
                    NodeScore += edge.NodeScore;
                    RatioScore += edge.RatioScore;
                    LcScore += edge.LcScore;
                    ImsScore += edge.ImsScore;
                }
            }
            Score = NodeScore + RatioScore + LcScore + ImsScore;
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

        private void UpdateEdges(IEnumerable<FeatureNode> nodes, List<FeatureNode> lNodes)
        {
            var maxWeight = double.NegativeInfinity;
            FeatureEdge edgeWithMaxWeight = null;
            
            foreach (var rnode in nodes)
            {
                foreach (var lnode in lNodes)
                {
                    if (rnode == lnode) continue;
                    var edge = new FeatureEdge(lnode, rnode, _scoringParams);
                    if (maxWeight < edge.Weight)
                    {
                        maxWeight = edge.Weight;
                        edgeWithMaxWeight = edge;
                    }
                }
            }
            if (edgeWithMaxWeight != null)
            {
                lNodes.Add(edgeWithMaxWeight.RNode);
                if(!ContainsKey(edgeWithMaxWeight.LNode)) this[edgeWithMaxWeight.LNode] = new List<FeatureEdge>(); 
                this[edgeWithMaxWeight.LNode].Add(edgeWithMaxWeight);
            } 
        }

        private List<FragmentFeatureNode> GetFragmentNodes(ImsDataCached imsData, Feature precursorFeature, Composition cutComposition, Ion precursorIon, GroupParameter parameter)
        {
            var ionTypes = _scoringParams.GetIonTypes(parameter);

            var nodes = new List<FragmentFeatureNode>();
            var suffixComposition = precursorIon.Composition -Composition.H2O - cutComposition;

            foreach (var ionType in ionTypes)
            {
                var composition = ionType.IsPrefixIon ? cutComposition : suffixComposition;
                //Console.WriteLine("FragFeatureGraph\t" + ionType.Name + "\t" + ionType.GetIon(composition).GetMz());
                var node = new FragmentFeatureNode(IsotopomerFeatures.GetFramentIsotopomerFeatures(imsData, composition, ionType, precursorFeature), ionType, parameter, _scoringParams);
                //if(node.Feature != null)
                nodes.Add(node);
            }
            return nodes;
        } 

    }
}
