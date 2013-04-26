using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.IMSScoring;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// ScoringGraph represents sequence variants with the same composition (e.g P*EPTIDE and PEP*TIDE)
    /// </summary>
    public class ScoringGraph
    {
        private const int DefaultMinPrecursorCharge = 1;
        private const int DefaultMaxPrecursorCharge = 4;

        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition _sequenceComposition;
        private readonly Dictionary<int,Ion> _precursorIon;

        private readonly ScoringGraphNode _rootNode;
        private readonly int _minPrecursorCharge;
        private readonly int _maxPrecursorCharge;

        private ImsDataCached _imsData;
        private ImsScorer _imsScorer;

        public ScoringGraph(AminoAcid[] aminoAcidSequence, Composition sequenceComposition, ScoringGraphNode rootNode)
            : this(aminoAcidSequence, sequenceComposition, rootNode, DefaultMinPrecursorCharge, DefaultMaxPrecursorCharge)
        {
        }

        public ScoringGraph(AminoAcid[] aminoAcidSequence, Composition sequenceComposition, ScoringGraphNode rootNode,
                            int minPrecursorCharge, int maxPrecursorCharge)
        {
            _aminoAcidSequence = aminoAcidSequence;
            _sequenceComposition = sequenceComposition;
            _rootNode = rootNode;
            _minPrecursorCharge = minPrecursorCharge;
            _maxPrecursorCharge = maxPrecursorCharge;

            _precursorIon = new Dictionary<int, Ion>();
            for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
            {
                _precursorIon[precursorCharge] = new Ion(_sequenceComposition, precursorCharge);
            } 
        }

        public void RegisterImsData(ImsDataCached imsData)
        {
            _imsData = imsData;
        }

        public void ComputeScores()
        {
            var compositions = new HashSet<Composition>();

            var curNodes = new HashSet<ScoringGraphNode> { _rootNode };
            int index = -1;
            while (curNodes.Any())
            {
                ++index;
                var newNodes = new HashSet<ScoringGraphNode>();
                foreach (var node in curNodes)
                {
                    compositions.Add(node.Composition);
                    foreach (var nextNode in node.GetNextNodes())
                    {
                        newNodes.Add(nextNode);
                    }
                }
                curNodes = newNodes;
            }
        }

        public IEnumerable<Composition> GetCompositions()
        {
            var compositions = new HashSet<Composition>();

            var curNodes = new HashSet<ScoringGraphNode> {_rootNode};
            while (curNodes.Any())
            {
                var newNodes = new HashSet<ScoringGraphNode>();
                foreach (var node in curNodes)
                {
                    compositions.Add(node.Composition);
                    foreach (var nextNode in node.GetNextNodes())
                    {
                        newNodes.Add(nextNode);
                    }
                }
                curNodes = newNodes;
            }

            return compositions;
        }

        private double GetCutScore(int index, Composition cutComposition, int precursorCharge)
        {
            Ion precursorIon = _precursorIon[precursorCharge];
            _imsData.GetPrecursorFeatures(precursorIon.GetMz());

            char nTermAA = _aminoAcidSequence[index - 1].Residue;
            char cTermAA = _aminoAcidSequence[index].Residue;
            //_imsScorer.GetCutScore(nTermAA, cTermAA, cutComposition, precursorIon);

            throw new System.NotImplementedException();
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
