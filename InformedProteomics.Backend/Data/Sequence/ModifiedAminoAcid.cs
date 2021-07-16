namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// A modified amino acid
    /// </summary>
    public class ModifiedAminoAcid : AminoAcid
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="aa"></param>
        /// <param name="modification"></param>
        public ModifiedAminoAcid(AminoAcid aa, Modification modification)
            : base(aa.Residue, aa.Name + "+" + modification.Name, aa.Composition + modification.Composition)
        {
            if (aa is not ModifiedAminoAcid modAa)
            {
                // aa is not modified
                Modification = modification;
            }
            else
            {
                // aa is already modified
                Modification = Modification.RegisterAndGetModification(modAa.Modification.Name + "+" + modification.Name, Composition);
            }
        }

        /// <summary>
        /// Modification applied to this amino acid
        /// </summary>
        public Modification Modification { get; }
    }
}
