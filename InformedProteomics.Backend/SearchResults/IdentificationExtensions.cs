using System;
using InformedProteomics.Backend.Data.Sequence;
using PSI_Interface.IdentData;

namespace InformedProteomics.Backend.SearchResults
{
    /// <summary>
    /// Extension functions for working with identification data.
    /// </summary>
    public static class IdentificationExtensions
    {
        /// <summary>
        /// Convert the sequence information from the external types to the internal types
        /// </summary>
        /// <param name="peptide"></param>
        /// <returns>Sequence object</returns>
        public static Sequence GetIpSequence(this SimpleMZIdentMLReader.PeptideRef peptide)
        {
            var aminoAcidSet = new AminoAcidSet();
            var sequence = new Sequence(peptide.Sequence, aminoAcidSet);
            var sequenceLength = sequence.Count;

            foreach (var mod in peptide.Mods)
            {
                // Mod.Key is the 1-based position of the mod in the peptide
                // However, N-terminal mods will have mod.Key = 0 and C-terminal mods will have mod.Key = peptide_length + 1

                // The following if statement assures that N- or C- terminal mods are associated with the first or last residue in the peptide, respectively

                int seqIndex;

                if (mod.Key < 1)
                    seqIndex = 0;
                else if (mod.Key >= sequenceLength)
                    seqIndex = sequenceLength - 1;
                else
                    seqIndex = mod.Key - 1;

                var mzidMod = mod.Value;

                var modification = Modification.Get(mzidMod.Tag, mzidMod.Mass);
                sequence[seqIndex] = new ModifiedAminoAcid(sequence[seqIndex], modification);
            }

            return new Sequence(sequence);
        }
    }
}
