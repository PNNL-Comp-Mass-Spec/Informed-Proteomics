using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class SequenceGraph
    {
        /// <summary>
        /// Create a graph representing the annotation. Annotation is reversed.
        /// </summary>
        /// <param name="aaSet">amino acid set</param>
        /// <param name="annotation">annotation (e.g. G.PEPTIDER.K or _.PEPTIDER._)</param>
        /// <returns></returns>
        public static SequenceGraph CreateGraph(AminoAcidSet aaSet, string annotation)
        {
            const char delimiter = (char)FastaDatabase.Delimiter;
            if (annotation == null || !Regex.IsMatch(annotation, @"^[A-Z" + delimiter + @"]\.[A-Z]+\.[A-Z" + delimiter + @"_]$")) return null;
            
            var nTerm = annotation[0] == FastaDatabase.Delimiter
                                  ? AminoAcid.ProteinNTerm
                                  : AminoAcid.PeptideNTerm;
            var cTerm = annotation[annotation.Length - 1] == FastaDatabase.Delimiter
                                  ? AminoAcid.ProteinCTerm
                                  : AminoAcid.PeptideCTerm;

            var sequence = annotation.Substring(2, annotation.Length - 4);
            return CreateGraph(aaSet, nTerm, sequence, cTerm);
        }

        /// <summary>
        /// Create a graph representing the sequence. Sequence is reversed.
        /// </summary>
        /// <param name="aaSet">amino acid set</param>
        /// <param name="sequence">sequence</param>
        /// <returns>sequence graph</returns>
        public static SequenceGraph CreatePeptideGraph(AminoAcidSet aaSet, string sequence)
        {
            var seqGraph = new SequenceGraph(aaSet, sequence.Length);
            seqGraph.AddAminoAcid(AminoAcid.PeptideCTerm.Residue);
            for (var i = sequence.Length - 1; i >= 0; i--)
            {
                if (seqGraph.AddAminoAcid(sequence[i]) == false) return null;
            }
            seqGraph.AddAminoAcid(AminoAcid.PeptideNTerm.Residue);
            return seqGraph;
        }

        /// <summary>
        /// Create a graph representing the sequence. Sequence is reversed.
        /// </summary>
        /// <param name="aaSet">amino acid set</param>
        /// <param name="nTerm">N-term amino acid</param>
        /// <param name="sequence">sequence</param>
        /// <param name="cTerm">C-term amino acid</param>
        /// <returns>sequence graph</returns>
        public static SequenceGraph CreateGraph(AminoAcidSet aaSet, AminoAcid nTerm, string sequence, AminoAcid cTerm)
        {
            var seqGraph = new SequenceGraph(aaSet, sequence.Length);
            seqGraph.AddAminoAcid(cTerm.Residue);
            for (var i = sequence.Length - 1; i >= 0; i--)
            {
                if (seqGraph.AddAminoAcid(sequence[i]) == false) return null;
            }
            seqGraph.AddAminoAcid(nTerm.Residue);
            return seqGraph;
        }

        public AminoAcidSet AminoAcidSet 
        {
            get { return _aminoAcidSet; } 
        }
        public ModificationParams ModificationParams 
        { 
            get { return _modificationParams; }
        }

        /// <summary>
        /// Gets the number of possible compositions of the current sequence 
        /// </summary>
        /// <returns>the number of possible compositions</returns>
        public int GetNumSequenceCompositions()
        {
            return _graph[_index].Length;
        }

        /// <summary>
        /// Gets the number of distinct compositions of the current sequence 
        /// </summary>
        /// <returns>the number of possible compositions</returns>
        public int GetNumDistinctSequenceCompositions()
        {
            var compositions = new HashSet<Composition.Composition>();
            for (var nodeIndex = 0; nodeIndex < _graph[_index].Length; nodeIndex++)
            {
                compositions.Add(GetComposition(_index, nodeIndex));
            }
            return compositions.Count;
        }
        
        public int GetNumSequenceCompositionsWithNTermCleavage(int numNTermCleavages)
        {
            return numNTermCleavages >= _index - 3 ? 0 : _graph[_index-1 - numNTermCleavages].Length;
        }

        /// <summary>
        /// Gets all possible compositions of the current sequence
        /// </summary>
        /// <returns>all possible compositions</returns>
        public Composition.Composition[] GetSequenceCompositions()
        {
            var numCompositions = _graph[_index].Length;
            var compositions = new Composition.Composition[numCompositions];
            for (var modIndex = 0; modIndex < numCompositions; modIndex++)
            {
                compositions[modIndex] = GetComposition(_index, modIndex);
            }
            return compositions;
        }

        public Composition.Composition[] GetSequenceCompositionsWithNTermCleavage(int numNTermCleavages)
        {
            if (numNTermCleavages >= _index - 3) return new Composition.Composition[0];

            var numCompositions = _graph[_index-1-numNTermCleavages].Length;
            var compositions = new Composition.Composition[numCompositions];
            for (var modIndex = 0; modIndex < numCompositions; modIndex++)
            {
                compositions[modIndex] = GetComposition(_index-1 - numNTermCleavages, modIndex);
            }
            return compositions;
        }
        
        /// <summary>
        /// Gets the composition of the sequence without variable modification.
        /// </summary>
        /// <returns></returns>
        public Composition.Composition GetUnmodifiedSequenceComposition()
        {
            return _prefixComposition[_index];
        }

        /// <summary>
        /// Gets the number of possible product compositions of the current sequence
        /// </summary>
        /// <returns>the number of possible product compositions</returns>
        public int GetNumFragmentCompositions()
        {
            var numNodes = 0;
            for (var index = 2; index < _graph.Length - 2; index++)
            {
                numNodes += _graph[index].Length;
            }
            return numNodes;
        }

        public IEnumerable<Composition.Composition> GetAllFragmentNodeCompositions()
        {
            for (var seqIndex = 2; seqIndex < _maxSeqIndex - 2; seqIndex++)
            {
                for (var modIndex = 0; modIndex < _graph[seqIndex].Length; modIndex++)
                {
                    yield return GetComposition(seqIndex, modIndex);
                }
            }
        }

        public IEnumerable<Composition.Composition> GetFragmentCompositions(int modIndex, int numNTermCleavages)
        {
            var seqIndex = _index - 1 - numNTermCleavages;

            var fragmentCompositions = new List<Composition.Composition>();

            var modIndexList = new HashSet<int> { modIndex };
            for (var si = seqIndex; si >= 2; si--)
            {
                var newModIndexList = new HashSet<int>();
                foreach (var mi in modIndexList)
                {
                    if (si < seqIndex) fragmentCompositions.Add(GetComposition(si, mi));
                    var node = _graph[si][mi];
                    foreach (var prevModIndex in node.GetPrevNodeIndices())
                    {
                        newModIndexList.Add(prevModIndex);
                    }
                }
                modIndexList = newModIndexList;
            }

            return fragmentCompositions;
        }

        public void SetSink(int modIndex, int numNTermCleavages)
        {
            _sinkSeqIndex = _index - 1 - numNTermCleavages;
            _sink = _graph[_sinkSeqIndex][modIndex];
            _sinkSequenceComposition = GetComposition(_sinkSeqIndex, modIndex);
            _sinkSequenceCompositionWithH2O = _sinkSequenceComposition + Composition.Composition.H2O;
            _sinkSequenceCompositionWithH2O.ComputeApproximateIsotopomerEnvelop();
            _compNodeComposition = new Composition.Composition[_maxSeqIndex,_modificationParams.NumModificationCombinations];
        }

        public Composition.Composition GetSinkSequenceCompositionWithH2O()
        {
            return _sinkSequenceCompositionWithH2O;
        }

        private Composition.Composition _sinkSequenceComposition;
        private Composition.Composition _sinkSequenceCompositionWithH2O;
        private Node _sink;
        private int _sinkSeqIndex;
        private Composition.Composition[,] _compNodeComposition;

        public double GetScore(int charge, IScorer scorer)
        {
            var precursorIon = new Ion(_sinkSequenceCompositionWithH2O, charge);

            // Get 
            var precursorIonScore = scorer.GetPrecursorIonScore(precursorIon);

            var nodeScore = new double?[_maxSeqIndex][];
            var maxScore = new double?[_maxSeqIndex][];
            for (var si = 0; si < _maxSeqIndex; si++)
            {
                nodeScore[si] = new double?[_graph[si].Length];
                maxScore[si] = new double?[_graph[si].Length];
                //if (si <= 1)
                //{
                //    for (var ni = 0; ni < nodeScore[si].Length; ni++) nodeScore[si][ni] = 0.0;
                //}
                //else
                //{
                //    for (var mi = 0; mi < nodeScore[si].Length; mi++)
                //    {
                //        var suffixComposition = GetComposition(si, mi);
                //        var prefixComposition = GetComplementaryComposition(si, mi);
                //        nodeScore[si][mi] = scorer.GetFragmentScore(prefixComposition, suffixComposition);
                //    }
                //}
            }
            maxScore[0][0] = 0;

            var fragmentScore =
                _sink.GetPrevNodeIndices()
                    .Max(prevNodeIndex => GetFragmentScore(_sinkSeqIndex - 1, prevNodeIndex, precursorIon, scorer, nodeScore, maxScore));
            return precursorIonScore + fragmentScore;
        }

        private double GetFragmentScore(int seqIndex, int modIndex, Ion precursorIon, IScorer scorer, double?[][] nodeScore, double?[][] maxScore)
        {
            var score = maxScore[seqIndex][modIndex];
            if (score != null) return (double)score;

            var node = _graph[seqIndex][modIndex];
            var curNodeScore = nodeScore[seqIndex][modIndex] ??
                (nodeScore[seqIndex][modIndex] = scorer.GetFragmentScore(GetComplementaryComposition(seqIndex, modIndex), GetComposition(seqIndex, modIndex)));

            var prevNodeScore = 0.0;
            if (node.GetPrevNodeIndices().Any())
            {
                prevNodeScore =
                    node.GetPrevNodeIndices()
                        .Max(prevNodeIndex => GetFragmentScore(seqIndex - 1, prevNodeIndex, precursorIon, scorer, nodeScore, maxScore));
            }
            maxScore[seqIndex][modIndex] = curNodeScore + prevNodeScore;
// ReSharper disable PossibleInvalidOperationException
            return (double)maxScore[seqIndex][modIndex];
// ReSharper restore PossibleInvalidOperationException
        }

        private Composition.Composition GetComposition(int seqIndex, int modIndex)
        {
            if (_nodeComposition[seqIndex, modIndex] == null)
            {
                var node = _graph[seqIndex][modIndex];
                _nodeComposition[seqIndex, modIndex] = _prefixComposition[seqIndex] +
                                  _modificationParams.GetModificationCombination(node.ModificationCombinationIndex)
                                                     .Composition;
            }
            return _nodeComposition[seqIndex, modIndex];
        }

        private Composition.Composition GetComplementaryComposition(int seqIndex, int modIndex)
        {
            if (_compNodeComposition[seqIndex, modIndex] == null)
            {
                var nodeComposition = GetComposition(seqIndex, modIndex);
                _compNodeComposition[seqIndex, modIndex] = _sinkSequenceComposition - nodeComposition;
            }
            return _compNodeComposition[seqIndex, modIndex];
        }

        private readonly int _maxSeqIndex;
        private readonly AminoAcidSet _aminoAcidSet;
        private readonly ModificationParams _modificationParams;

        private int _index;
        private readonly Node[][] _graph;
        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition.Composition[,] _nodeComposition;

        private readonly Composition.Composition[] _prefixComposition;

        /// <summary>
        /// Initialize and set the maximum sequence length
        /// </summary>
        /// <param name="aminoAcidSet"></param>
        /// <param name="maxSequenceLength"></param>
        private SequenceGraph(AminoAcidSet aminoAcidSet, int maxSequenceLength)
        {
            _aminoAcidSet = aminoAcidSet;
            _modificationParams = aminoAcidSet.GetModificationParams();

            _maxSeqIndex = maxSequenceLength + 3;  // init + C-term + sequence length + N-term

            _index = 0;

            _aminoAcidSequence = new AminoAcid[_maxSeqIndex];
            _aminoAcidSequence[0] = AminoAcid.Empty;

            _prefixComposition = new Composition.Composition[_maxSeqIndex];
            _prefixComposition[0] = Composition.Composition.Zero;

            _graph = new Node[_maxSeqIndex][];
            _graph[0] = new[] { new Node(0) };

            _nodeComposition = new Composition.Composition[_maxSeqIndex,_modificationParams.NumModificationCombinations];
        }

        private bool AddAminoAcid(char residue)
        {
            return PutAminoAcid(_index, residue);
        }

        /// <summary>
        /// Add an amino acid residue to this generator.
        /// </summary>
        /// <param name="index">index to add the amino acid. 0 is N-term. 1 is the first amino acid.</param>
        /// <param name="residue">amino acid residue to add.</param>
        /// <returns>true if residue is a valid amino acid; false otherwise.</returns>
        private bool PutAminoAcid(int index, char residue)
        {
            _index = index + 1;

            var aminoAcid = AminoAcidSet.GetAminoAcid(residue);
            if (aminoAcid == null) // residue is not valid
            {
                return false;
            }

            _aminoAcidSequence[_index] = aminoAcid;
            _prefixComposition[_index] = _prefixComposition[_index - 1] + aminoAcid.Composition;

            // TODO: Location-specific modifications (e.g. N-term Acetyl) are not supported yet
            var modIndices = AminoAcidSet.GetModificationIndices(residue);

            if (!modIndices.Any())  // No modification
            {
                _graph[_index] = new Node[_graph[_index - 1].Length];
                for (var i = 0; i < _graph[_index - 1].Length; i++)
                {
                    _graph[_index][i] = new Node(_graph[_index - 1][i].ModificationCombinationIndex, i);
                }
            }
            else
            {
                var modCombIndexToNodeMap = new Dictionary<int, Node>();
                for (var i = 0; i < _graph[_index - 1].Length; i++)
                {
                    var prevNodeIndex = i;
                    var prevNode = _graph[_index - 1][i];
                    var prevModCombIndex = prevNode.ModificationCombinationIndex;
                    Node newNode;
                    // unmodified edge
                    if (modCombIndexToNodeMap.TryGetValue(prevModCombIndex, out newNode))
                    {
                        newNode.AddPrevNodeIndex(prevNodeIndex);
                    }
                    else
                    {
                        modCombIndexToNodeMap.Add(prevModCombIndex, new Node(prevModCombIndex, prevNodeIndex));
                    }

                    // modified edge
                    foreach (var modIndex in modIndices)
                    {
                        var modCombIndex = ModificationParams.GetModificationCombinationIndex(
                                                    prevNode.ModificationCombinationIndex, modIndex);
                        if (modCombIndex < 0)   // too many modifications
                            continue;
                        if (modCombIndexToNodeMap.TryGetValue(modCombIndex, out newNode))
                        {
                            newNode.AddPrevNodeIndex(prevNodeIndex);
                        }
                        else
                        {
                            modCombIndexToNodeMap.Add(modCombIndex, new Node(modCombIndex, prevNodeIndex));
                        }
                    }
                    _graph[_index] = modCombIndexToNodeMap.Values.ToArray();
                }
            }

            return true;
        }

        //public IEnumerable<ScoringGraph> GetScoringGraphs()
        //{
        //    var numCompositions = _graph[_index].Length;
        //    var scoringGraphs = new ScoringGraph[numCompositions];
        //    for (var i = 0; i < numCompositions; i++)
        //    {
        //        scoringGraphs[i] = GetScoringGraph(i);
        //    }
        //    return scoringGraphs;
        //}

        //public ScoringGraph GetScoringGraph(int sequenceIndex)
        //{
        //    // backtracking
        //    ScoringGraphNode rootNode = null;
        //    var nextNodeMap = new Dictionary<int, List<ScoringGraphNode>>
        //        {
        //            {sequenceIndex, new List<ScoringGraphNode>()}
        //        };
        //    for (int modIndex = _index; modIndex >= 0; --modIndex )
        //    {
        //        var newNextNodeMap = new Dictionary<int, List<ScoringGraphNode>>();

        //        foreach (var entry in nextNodeMap)
        //        {
        //            int modIndex = entry.Key;
        //            var curNode = _graph[modIndex][modIndex];
        //            var composition = GetComposition(modIndex, modIndex);
        //            var scoringGraphNode = new ScoringGraphNode(composition, modIndex);
        //            scoringGraphNode.AddNextNodes(nextNodes: entry.Value);

        //            foreach (var prevNodeIndex in curNode.GetPrevNodeIndices())
        //            {
        //                List<ScoringGraphNode> nextNodes;
        //                if (newNextNodeMap.TryGetValue(prevNodeIndex, out nextNodes))
        //                {
        //                    nextNodes.Add(scoringGraphNode);
        //                }
        //                else
        //                {
        //                    newNextNodeMap.Add(prevNodeIndex, new List<ScoringGraphNode> { scoringGraphNode });
        //                }
        //            }

        //            if (!curNode.GetPrevNodeIndices().Any())
        //            {
        //                rootNode = scoringGraphNode;
        //            }
        //        }
                
        //        nextNodeMap = newNextNodeMap;
        //    }
            
        //    var scoringGraph = new ScoringGraph(_aminoAcidSequence, GetComposition(_index, sequenceIndex), rootNode);

        //    return scoringGraph;
        //}

    }

    internal class Node
    {
        internal Node(int modificationCombinationIndex)
        {
            ModificationCombinationIndex = modificationCombinationIndex;
            _prevNodeIndices = new int[2];
            _count = 0;
        }

        internal Node(int modificationCombinationIndex, int prevNodeIndex)
            : this(modificationCombinationIndex)
        {
            AddPrevNodeIndex(prevNodeIndex);
        }

        public int ModificationCombinationIndex { get; private set; }
        
        private int[] _prevNodeIndices;
        private int _count;

        internal bool AddPrevNodeIndex(int prevNodeIndex)
        {
            for (var i = 0; i < _count; i++)
            {
                if (_prevNodeIndices[i] == prevNodeIndex) return false;
            }
            if(_count >= _prevNodeIndices.Length)
            {
                Array.Resize(ref _prevNodeIndices, _prevNodeIndices.Length * 2); 
            }
            _prevNodeIndices[_count++] = prevNodeIndex; 

            return true;
        }

        public IEnumerable<int> GetPrevNodeIndices()
        {
            //return _prevNodeIndices;
            for (var i = 0; i < _count; i++)
            {
                yield return _prevNodeIndices[i];
            }
        }
    }

}
