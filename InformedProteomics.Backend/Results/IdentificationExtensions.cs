namespace InformedProteomics.Backend.Results
{
    using System;

    using InformedProteomics.Backend.Data.Sequence;

    using PSI_Interface.IdentData;

    /// <summary>
    /// Extension functions for working with identification data.
    /// </summary>
    public static class IdentificationExtensions
    {
        /// <summary>
        /// Convert the sequence information from the external types to the internal types
        /// </summary>
        /// <param name="peptide"></param>
        /// <returns></returns>
        public static Sequence GetIpSequence(this SimpleMZIdentMLReader.PeptideRef peptide)
        {
            var aminoAcidSet = new AminoAcidSet();
            var sequence = new Sequence(peptide.Sequence, aminoAcidSet);

            foreach (var mod in peptide.Mods)
            {
                var seqIndex = Math.Max(0, mod.Key - 1);
                var mzidMod = mod.Value;

                var modification = Modification.Get(mzidMod.Tag, mzidMod.Mass);
                sequence[seqIndex] = new ModifiedAminoAcid(sequence[seqIndex], modification);
            }

            return new Sequence(sequence);
        }
    }
}
