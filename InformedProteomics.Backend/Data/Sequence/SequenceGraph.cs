using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

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
        private readonly int _maxNumDynModsPerSequence;
        private readonly int _maxNumModificationCombinations;
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
            int index = 0;
            PutAminoAcid(index, AminoAcid.PeptideNTerm.Residue);
            foreach (char aaResidue in pepSequence)
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

            _maxNumDynModsPerSequence = _modificationParams.MaxNumDynModsPerSequence;
            _maxNumModificationCombinations = _modificationParams.NumModificationCombinations;

            _maxIndex = maxSequenceLength + 3;  // init + N-term + sequence length + C-term

            _index = 0;

            _aminoAcidSequence = new AminoAcid[_maxIndex];
            _aminoAcidSequence[0] = AminoAcid.Empty;

            _prefixComposition = new Composition[_maxIndex];
            _prefixComposition[0] = Composition.Zero;

            _graph = new Node[_maxIndex][];
            _graph[0] = new[] { new Node(0) };
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

            AminoAcid aminoAcid = AminoAcidSet.GetAminoAcid(residue);
            if (aminoAcid == null) // residue is not valid
                return false;

            _aminoAcidSequence[_index] = aminoAcid;
            _prefixComposition[_index] = _prefixComposition[_index - 1] + aminoAcid.Composition;

            // TODO: Support location-specific modifications
            var modIndices = AminoAcidSet.GetModificationIndices(residue); 

            if (!modIndices.Any())  // No modification
            {
                _graph[_index] = new Node[_graph[_index - 1].Length];
                for (int i = 0; i < _graph[_index - 1].Length; i++)
                {
                    _graph[_index][i] = new Node(_graph[_index - 1][i].ModificationCombinationIndex, _graph[_index][i]);
                }
            }
            else
            {
                var modCombIndexToNodeMap = new Dictionary<int, Node>();
                for (int i = 0; i < _graph[_index - 1].Length; i++)
                {
                    Node prevNode = _graph[_index - 1][i];
                    int prevModCombIndex = prevNode.ModificationCombinationIndex;
                    Node newNode;
                    // unmodified
                    if(modCombIndexToNodeMap.TryGetValue(prevModCombIndex, out newNode))
                    {
                        newNode.AddPrevNode(prevNode);
                    }
                    else
                    {
                        modCombIndexToNodeMap.Add(prevModCombIndex, new Node(prevModCombIndex, prevNode));
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
                            newNode.AddPrevNode(prevNode);
                        }
                        else
                        {
                            modCombIndexToNodeMap.Add(modCombIndex, new Node(modCombIndex, prevNode));
                        }
                    }
                    _graph[_index] = modCombIndexToNodeMap.Values.ToArray();
                }
            }

            return true;
        }

        /// <summary>
        /// Gets all possible compositions of the current sequence
        /// </summary>
        /// <returns></returns>
        public Composition[] GetSequenceCompositions()
        {
            Composition unmodComp = GetUnmodifiedSequenceComposition();

            var compositions = new List<Composition>();
            compositions.AddRange(_graph[_index].Select(
                node => unmodComp + _modificationParams.GetModificationCombination(node.ModificationCombinationIndex).Composition));
            return compositions.ToArray();
        }

        /// <summary>
        /// Gets the composition of the sequence without variable modification.
        /// </summary>
        /// <returns></returns>
        public Composition GetUnmodifiedSequenceComposition()
        {
            return _prefixComposition[_index];
        }

        /// <summary>
        /// Gets possible compositions of a prefix fragment.
        /// index=0 means 1st cleavage (e.g. P/EPTIDE) 
        /// </summary>
        /// <param name="sequenceIndex">Index of the sequence. 0 for sequence with no variable modification.</param>
        /// <returns></returns>
        public Composition[][] GetFragmentCompositions(int sequenceIndex)
        {
           // // backtracking
           // var fragmentCompositions = new Composition[_index - 2][];

           // var activeNodes = new List<Node> {_graph[_index][sequenceIndex]};
           // for (int index = _index - 1; index >= 2; --index)
           // {
           //     var newActiveNodes = new List<Node>();
           //     foreach (var node in activeNodes)
           //     {
                    
           //     }
           // }
           //return fragmentCompositions;
            throw new System.NotImplementedException();
        }
    }

    internal class Node
    {
        internal Node(int modificationCombinationIndex)
        {
            ModificationCombinationIndex = modificationCombinationIndex;
            _prevNodes = new List<Node>();
        }

        internal Node(int modificationCombinationIndex, Node prevNode): this(modificationCombinationIndex)
        {
            _prevNodes.Add(prevNode);
        }

        internal int ModificationCombinationIndex { get; private set; }

        private readonly IList<Node> _prevNodes;

        internal void AddPrevNode(Node prevNode)
        {
            _prevNodes.Add(prevNode);
        }

        internal IEnumerable<Node> GetPrevNodes() { return _prevNodes; }
    }
}
