using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Sequence
{
    class SequenceGenerator
    {
        private const int MaxSequenceLengthForBottomUp = 50;

        public AminoAcidSet AminoAcidSet { get; private set; }
        private readonly int _maxIndex;
        private readonly int _maxNumDynModsPerPeptide;

        private int _index;
        private Composition[][] seqComposition;

        public SequenceGenerator(AminoAcidSet aminoAcidSet): this(aminoAcidSet, MaxSequenceLengthForBottomUp)
        {
        }

        public SequenceGenerator(AminoAcidSet aminoAcidSet, int maxSequenceLength)
        {
            AminoAcidSet = aminoAcidSet;
            _maxNumDynModsPerPeptide = aminoAcidSet.MaxNumDynModsPerPeptide;
            _maxIndex = maxSequenceLength + 3;  // init + N-term + sequence length + C-term

            _index = 0;
            seqComposition = new Composition[_maxIndex][];
            seqComposition[0] = new[] { Composition.Zero };
        }

        /// <summary>
        /// Add an amino acid residue to this generator.
        /// </summary>
        /// <param name="residue"></param>
        /// <returns>true if residue is a valid amino acid; false otherwise.</returns>
        public bool AddAminoAcid(char residue)
        {
            ++_index;
            AminoAcid[] aminoAcids = AminoAcidSet.GetAminoAcids(residue);
            if (aminoAcids == null) // residue is not valid
                return false;

    
            return true;
        }

        /// <summary>
        /// Gets possible compositions of a prefix fragment.
        /// index=0 means N-terminus. index=i (0 < i 
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public Composition[] GetFragmentCompositions(int index)
        {
            throw new System.NotImplementedException();
        }

    }
}
