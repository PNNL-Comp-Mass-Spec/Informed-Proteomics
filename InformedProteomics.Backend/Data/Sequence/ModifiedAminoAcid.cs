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
            var modAa = aa as ModifiedAminoAcid;
            if (modAa == null) _modification = modification; // aa is not modified
            else    // aa is already modified
            {
                _modification = Modification.RegisterAndGetModification(modAa.Modification.Name + "+" + modification.Name, Composition);
            }
        }

        private readonly Modification _modification;

        /// <summary>
        /// Modification applied to this amino acid
        /// </summary>
        public Modification Modification
        {
            get { return _modification; }
        }
    }
}
