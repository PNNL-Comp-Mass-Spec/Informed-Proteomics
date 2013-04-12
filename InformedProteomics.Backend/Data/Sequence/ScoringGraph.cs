using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.IMSScoring;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class ScoringGraph
    {
        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition _sequenceComposition;
        private readonly ScoringGraphNode _rootNode;
        private ImsScorer _imsScorer;

        public ScoringGraph(AminoAcid[] aminoAcidSequence, Composition sequenceComposition, ScoringGraphNode rootNode)
        {
            _aminoAcidSequence = aminoAcidSequence;
            _sequenceComposition = sequenceComposition;
            _rootNode = rootNode;
        }

        public void RegisterScorer(ImsScorer imsScorer)
        {
            _imsScorer = imsScorer;
        }

        public void ComputeScores()
        {
            throw new System.NotImplementedException();
        }

        public IEnumerable<Composition> GetCompositions()
        {
            return GetCompositions(_rootNode);
        }

        public ISet<Composition> GetCompositions(ScoringGraphNode node)
        {
            var compositions = new HashSet<Composition> {node.Composition};
            foreach (var nextNode in node.GetNextNodes())
            {
                compositions.Add(nextNode.Composition);
            }

            return compositions;
        }
    }

    public class ScoringGraphNode
    {
        private readonly Composition _composition;
        private readonly List<ScoringGraphNode> _nextNodes;

        public ScoringGraphNode(Composition composition)
        {
            _composition = composition;
            _nextNodes = new List<ScoringGraphNode>();
        }

        public Composition Composition
        {
            get { return _composition; }
        }

        public void AddNextNode(ScoringGraphNode nextNode)
        {
            _nextNodes.Add(nextNode);
        }

        public void AddNextNodes(IEnumerable<ScoringGraphNode> nextNodes)
        {
            _nextNodes.AddRange(nextNodes);
        }

        public IEnumerable<ScoringGraphNode> GetNextNodes()
        {
            return _nextNodes;
        }

        public double Score { get; internal set; }
    }
}
