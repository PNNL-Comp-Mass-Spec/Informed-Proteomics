using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class SequenceGraph
    {
        private const int MaxPeptideLength = 50;

        public AminoAcidSet AminoAcidSet 
        {
            get { return _aminoAcidSet; } 
        }
        public ModificationParams ModificationParams 
        { 
            get { return _modificationParams; }
        }

        private readonly int _maxIndex;
        private readonly AminoAcidSet _aminoAcidSet;
        private readonly ModificationParams _modificationParams;

        private int _index;
        private readonly Node[][] _graph;
        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition[] _prefixComposition;

        /// <summary>
        /// Initialize. Set the maximum sequence length as MaxPeptideLength
        /// </summary>
        /// <param name="aminoAcidSet"></param>
        public SequenceGraph(AminoAcidSet aminoAcidSet)
            : this(aminoAcidSet, MaxPeptideLength)
        {
        }

        /// <summary>
        /// Initialize using a peptide sequence
        /// </summary>
        /// <param name="aminoAcidSet">amino acid set</param>
        /// <param name="pepSequence">peptide sequence</param>
        public SequenceGraph(AminoAcidSet aminoAcidSet, string pepSequence)
            : this(aminoAcidSet, pepSequence.Length)
        {
            var index = 0;
            PutAminoAcid(index, AminoAcid.PeptideNTerm.Residue);
            foreach (var aaResidue in pepSequence)
            {
                ++index;
                PutAminoAcid(index, aaResidue);
            }
            PutAminoAcid(++index, AminoAcid.PeptideCTerm.Residue);
        }

        /// <summary>
        /// Initialize and set the maximum sequence length
        /// </summary>
        /// <param name="aminoAcidSet"></param>
        /// <param name="maxSequenceLength"></param>
        public SequenceGraph(AminoAcidSet aminoAcidSet, int maxSequenceLength)
        {
            _aminoAcidSet = aminoAcidSet;
            _modificationParams = aminoAcidSet.GetModificationParams();

            _maxIndex = maxSequenceLength + 3;  // init + N-term + sequence length + C-term

            _index = 0;

            _aminoAcidSequence = new AminoAcid[_maxIndex];
            _aminoAcidSequence[0] = AminoAcid.Empty;

            _prefixComposition = new Composition[_maxIndex];
            _prefixComposition[0] = Composition.Zero;

            _graph = new Node[_maxIndex][];
            _graph[0] = new[] { new Node(0) };
        }

        public bool PutPeptideSequence(string sequence)
        {
            return PutSequence(sequence, false);
        }

        public bool PutSequence(string sequence, bool isProtein)
        {
            var index = 0;
            PutAminoAcid(index, isProtein ? AminoAcid.ProteinNTerm.Residue : AminoAcid.PeptideNTerm.Residue);

            foreach (var aaResidue in sequence)
            {
                ++index;
                if (PutAminoAcid(index, aaResidue) == false) return false;
            }
            return true;
        }

        public bool AddAminoAcid(char residue)
        {
            return PutAminoAcid(_index, residue);
        }

        /// <summary>
        /// Add an amino acid residue to this generator.
        /// </summary>
        /// <param name="index">index to add the amino acid. 0 is N-term. 1 is the first amino acid.</param>
        /// <param name="residue">amino acid residue to add.</param>
        /// <returns>true if residue is a valid amino acid; false otherwise.</returns>
        public bool PutAminoAcid(int index, char residue)
        {
            _index = index + 1;

            var aminoAcid = AminoAcidSet.GetAminoAcid(residue);
            if (aminoAcid == null) // residue is not valid
            {
                return false;
            }

            _aminoAcidSequence[_index] = aminoAcid;
            _prefixComposition[_index] = _prefixComposition[_index - 1] + aminoAcid.Composition;

            // TODO: Support location-specific modifications
            var modIndices = AminoAcidSet.GetModificationIndices(residue); 

            if (!modIndices.Any())  // No modification
            {
                _graph[_index] = new Node[_graph[_index - 1].Length];
                for (int i = 0; i < _graph[_index - 1].Length; i++)
                {
                    _graph[_index][i] = new Node(_graph[_index - 1][i].ModificationCombinationIndex, i);
                }
            }
            else
            {
                var modCombIndexToNodeMap = new Dictionary<int, Node>();
                for (int i = 0; i < _graph[_index - 1].Length; i++)
                {
                    int prevNodeIndex = i;
                    Node prevNode = _graph[_index - 1][i];
                    int prevModCombIndex = prevNode.ModificationCombinationIndex;
                    Node newNode;
                    // unmodified
                    if(modCombIndexToNodeMap.TryGetValue(prevModCombIndex, out newNode))
                    {
                        newNode.AddPrevNodeIndex(prevNodeIndex);
                    }
                    else
                    {
                        modCombIndexToNodeMap.Add(prevModCombIndex, new Node(prevModCombIndex, prevNodeIndex));
                    }
                    // modified
                    foreach(var modIndex in modIndices)
                    {
                        int modCombIndex = ModificationParams.GetModificationCombinationIndex(
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

        /// <summary>
        /// Gets the number of possible compositions of the current sequence 
        /// </summary>
        /// <returns>the number of possible compositions</returns>
        public int GetNumCompositions()
        {
            return _graph[_index].Length;
        }

        /// <summary>
        /// Gets all possible compositions of the current sequence
        /// </summary>
        /// <returns>all possible compositions</returns>
        public Composition[] GetSequenceCompositions()
        {
            var numCompositions = _graph[_index].Length;
            var compositions = new Composition[numCompositions];
            for (var nodeIndex = 0; nodeIndex < numCompositions; nodeIndex++)
            {
                compositions[nodeIndex] = GetComposition(_index, nodeIndex);
            }
            return compositions;
        }

        /// <summary>
        /// Gets the composition of the sequence without variable modification.
        /// </summary>
        /// <returns></returns>
        public Composition GetUnmodifiedSequenceComposition()
        {
            return _prefixComposition[_index];
        }

        public IEnumerable<ScoringGraph> GetScoringGraphs()
        {
            int numCompositions = _graph[_index].Length;
            var scoringGraphs = new ScoringGraph[numCompositions];
            for (int i = 0; i < numCompositions; i++)
            {
                scoringGraphs[i] = GetScoringGraph(i);
            }
            return scoringGraphs;
        }

        public ScoringGraph GetScoringGraph(int sequenceIndex)
        {
            // backtracking
            ScoringGraphNode rootNode = null;
            var nextNodeMap = new Dictionary<int, List<ScoringGraphNode>>
                {
                    {sequenceIndex, new List<ScoringGraphNode>()}
                };
            for (int seqIndex = _index; seqIndex >= 0; --seqIndex )
            {
                var newNextNodeMap = new Dictionary<int, List<ScoringGraphNode>>();

                foreach (var entry in nextNodeMap)
                {
                    int nodeIndex = entry.Key;
                    var curNode = _graph[seqIndex][nodeIndex];
                    var composition = GetComposition(seqIndex, nodeIndex);
                    var scoringGraphNode = new ScoringGraphNode(composition, seqIndex);
                    scoringGraphNode.AddNextNodes(nextNodes: entry.Value);

                    foreach (var prevNodeIndex in curNode.GetPrevIndices())
                    {
                        List<ScoringGraphNode> nextNodes;
                        if (newNextNodeMap.TryGetValue(prevNodeIndex, out nextNodes))
                        {
                            nextNodes.Add(scoringGraphNode);
                        }
                        else
                        {
                            newNextNodeMap.Add(prevNodeIndex, new List<ScoringGraphNode> { scoringGraphNode });
                        }
                    }

                    if (!curNode.GetPrevIndices().Any())
                    {
                        rootNode = scoringGraphNode;
                    }
                }
                
                nextNodeMap = newNextNodeMap;
            }
            
            var scoringGraph = new ScoringGraph(_aminoAcidSequence, GetComposition(_index, sequenceIndex), rootNode);

            return scoringGraph;
        }

        private Composition GetComposition(int seqIndex, int nodeIndex)
        {
            Node node = _graph[seqIndex][nodeIndex];
            var composition = _prefixComposition[seqIndex] +
                              _modificationParams.GetModificationCombination(node.ModificationCombinationIndex)
                                                 .Composition;
            return composition;
        }
    }

    internal class Node
    {
        internal Node(int modificationCombinationIndex)
        {
            ModificationCombinationIndex = modificationCombinationIndex;
            _prevNodeIndices = new HashSet<int>();
        }

        internal Node(int modificationCombinationIndex, int prevNodeIndex)
            : this(modificationCombinationIndex)
        {
            _prevNodeIndices.Add(prevNodeIndex);
        }

        public int ModificationCombinationIndex { get; private set; }

        private readonly HashSet<int> _prevNodeIndices;

        internal bool AddPrevNodeIndex(int prevNodeIndex)
        {
            return _prevNodeIndices.Add(prevNodeIndex);
        }

        public IEnumerable<int> GetPrevIndices() { return _prevNodeIndices; }
    }

}
