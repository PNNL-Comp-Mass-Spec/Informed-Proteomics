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
        private readonly ScoringGraphNode[] _nodes;

        private readonly int _minPrecursorCharge;
        private readonly int _maxPrecursorCharge;

        private ImsDataCached _imsData;
        private ImsScorerFactory _imsScorerFactory;

        internal ScoringGraph(AminoAcid[] aminoAcidSequence, Composition sequenceComposition, ScoringGraphNode rootNode)
            : this(aminoAcidSequence, sequenceComposition, rootNode, DefaultMinPrecursorCharge, DefaultMaxPrecursorCharge)
        {
        }

        internal ScoringGraph(AminoAcid[] aminoAcidSequence, Composition sequenceComposition, ScoringGraphNode rootNode,
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

            // store all nodes in an array
            var nodes = new HashSet<ScoringGraphNode>();

            var curNodes = new HashSet<ScoringGraphNode> { _rootNode };
            while (curNodes.Any())
            {
                var newNodes = new HashSet<ScoringGraphNode>();
                foreach (var node in curNodes)
                {
                    if (nodes.Add(node))    // if node is new
                    {
                        foreach (var nextNode in node.GetNextNodes())
                        {
                            newNodes.Add(nextNode);
                        }
                    }
                }
                curNodes = newNodes;
            }

            _nodes = nodes.ToArray();
        }

        public void RegisterImsData(ImsDataCached imsData, ImsScorerFactory imsScorerFactory)
        {
            _imsData = imsData;
            _imsScorerFactory = imsScorerFactory;
        }

        // Feature, precursorScore, totalScore
        public Tuple<Feature, double, double> GetBestFeatureAndScore(int precursorCharge)
        {
            var precursorIon = new Ion(_sequenceComposition + Composition.H2O, precursorCharge);
            var imsScorer = _imsScorerFactory.GetImsScorer(_imsData, precursorIon);

            var precursorFeatureSet = _imsData.GetPrecursorFeatures(precursorIon.GetMz());
            var bestPrecursorScore = double.NegativeInfinity;
            var bestScore = double.NegativeInfinity;
            Feature bestFeature = null;
            foreach (var precursorFeature in precursorFeatureSet)
            {
                var precursorScore = imsScorer.GetPrecursorScore(precursorFeature);
                var productScore = GetProductIonScore(imsScorer, precursorFeature);

                var curFeatureScore = precursorScore + productScore;
                if (curFeatureScore > bestScore)
                {
                    bestPrecursorScore = precursorScore;
                    bestScore = curFeatureScore;
                    bestFeature = precursorFeature;
                }
            }

            return new Tuple<Feature, double, double>(bestFeature, bestPrecursorScore, bestScore);
        }

        public Ion GetPrecursorIon(int charge)
        {
            return new Ion(_sequenceComposition+Composition.H2O, charge);
        }

        public IEnumerable<Composition> GetCompositions()
        {
            return _nodes.Select(node => node.Composition).ToArray();
        }

        public string GetUnmodSequence()
        {
            var pepSequence = new StringBuilder();
            foreach (var aa in _aminoAcidSequence)
            {
                pepSequence.Append(aa.Residue);
            }
            return pepSequence.ToString();
        }

        private double GetProductIonScore(ImsScorer imsScorer, Feature precursorFeature)
        {
            return GetProductIonScore(_rootNode, imsScorer, precursorFeature);
        }

        private double GetProductIonScore(ScoringGraphNode node, ImsScorer imsScorer, Feature precursorFeature)
        {
            //Console.Write("Index: " + node.Index);
            double cutScore = 0;
            if (node.Index > 1 && node.Index <= _aminoAcidSequence.Length-3)
            {
                var nTermAA = _aminoAcidSequence[node.Index].Residue;
                var cTermAA = _aminoAcidSequence[node.Index + 1].Residue;
                cutScore = imsScorer.GetCutScore(nTermAA, cTermAA, node.Composition, precursorFeature);
                //Console.WriteLine(" " + _aminoAcidSequence[node.Index].Residue + " " + node.Composition + " " + _aminoAcidSequence[node.Index + 1].Residue + " " + cutScore);
            }
            //Console.WriteLine();

            if (node.Index > _aminoAcidSequence.Length - 3)
                return 0;
            var nextNodeScore = node.GetNextNodes().DefaultIfEmpty().Max(nextNode => GetProductIonScore(nextNode, imsScorer, precursorFeature));
            return cutScore + nextNodeScore;
        }
    }

    public class ScoringGraphNode
    {
        private readonly Composition _composition;
        private readonly int _index;
        private readonly List<ScoringGraphNode> _nextNodes;

        public ScoringGraphNode(Composition composition, int index)
        {
            _composition = composition;
            _index = index;
            _nextNodes = new List<ScoringGraphNode>();
        }

        public Composition Composition
        {
            get { return _composition; }
        }

        public int Index
        {
            get { return _index; }
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
