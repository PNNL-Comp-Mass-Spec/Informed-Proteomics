using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class SequenceGraph
    {
        private const int MaxPeptideLength = 50;

        public AminoAcidSet AminoAcidSet { get; private set; }

        private readonly int _maxIndex;
        private readonly int _maxNumDynModsPerPeptide;

        private int _index;
        private Node[][] graph;

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
            AminoAcidSet = aminoAcidSet;
            _maxNumDynModsPerPeptide = aminoAcidSet.MaxNumDynModsPerPeptide;
            _maxIndex = maxSequenceLength + 3;  // init + N-term + sequence length + C-term

            _index = 0;

            graph = new Node[_maxIndex][];
            graph[0] = new[] { new Node(Composition.Zero) };
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

            AminoAcid[] aminoAcids = AminoAcidSet.GetAminoAcids(residue);
            if (aminoAcids == null) // residue is not valid
                return false;

            if (aminoAcids.Length == 1) // unmodified
            {
                AminoAcid aa = aminoAcids[0];
                graph[_index] = new Node[graph[_index - 1].Length];
                for (int i = 0; i < graph.Length; i++)
                {
                    graph[_index][i] = new Node(graph[_index - 1][i], aa.Composition);
                }
            }
            else
            {
                //var compositionSet = new HashSet<Composition>();
                //var newCompositionList = new List<Composition>();
                //foreach (var prevSeqComposition in seqComposition[_index - 1])
                //{
                //    newCompositionList.AddRange(aminoAcids.Select(aa => prevSeqComposition + aa.Composition).Where(compositionSet.Add));
                //}
                //seqComposition[_index] = newCompositionList.ToArray();
            }

            return true;
        }

        /// <summary>
        /// Gets all possible compositions of the current sequence
        /// </summary>
        /// <returns></returns>
        public Composition[] GetSequenceCompositions()
        {
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Gets the composition of the sequence without variable modification.
        /// </summary>
        /// <returns></returns>
        public Composition GetUnmodifiedSequenceComposition()
        {
            throw new System.NotImplementedException();
        }

        /// <summary>
        /// Gets possible compositions of a prefix fragment.
        /// index=0 means N-terminus. index=i (0 < i 
        /// </summary>
        /// <param name="sequenceIndex">Index of the sequence. 0 for sequence with no variable modification.</param>
        /// <returns></returns>
        public Composition[][] GetFragmentCompositions(int sequenceIndex)
        {
            throw new System.NotImplementedException();
        }


    }

    internal class Node
    {
        internal Node(Composition composition)
        {
            Composition = composition;
            _prevNodes = new List<Node>();
        }

        internal Node(Node prevNode, Composition composition)
            : this(prevNode.Composition + composition)
        {
            _prevNodes.Add(prevNode);
        }

        internal Composition Composition { get; private set; }

        private readonly IList<Node> _prevNodes;

        internal void AddPrevNode(Node prevNode)
        {
            _prevNodes.Add(prevNode);
        }

        internal IEnumerable<Node> GetPrevNodes() { return _prevNodes; }
    }
}
