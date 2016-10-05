namespace InformedProteomics.Backend.Results
{
    using InformedProteomics.Backend.Data.Sequence;

    using PSI_Interface.IdentData;

    public static class IdentificationExtensions
    {
        public static Sequence GetIpSequence(this SimpleMZIdentMLReader.PeptideRef peptide)
        {
            var aminoAcidSet = new AminoAcidSet();
            var sequence = new Sequence(peptide.Sequence, aminoAcidSet);

            foreach (var mod in peptide.Mods)
            {
                var seqIndex = mod.Key - 1;
                var mzidMod = mod.Value;

                var modification = Modification.Get(mzidMod.Tag);
                sequence[seqIndex] = new ModifiedAminoAcid(sequence[seqIndex], modification);
            }

            return new Sequence(sequence);
        }
    }
}
