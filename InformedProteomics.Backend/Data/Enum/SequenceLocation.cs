namespace InformedProteomics.Backend.Data.Enum
{
    /// <summary>
    /// Allowed locations for a modification within a sequence
    /// </summary>
    public enum SequenceLocation
    {
        /// <summary>
        /// Modification is allowed at any position in the sequence
        /// </summary>
        Everywhere,

        /// <summary>
        /// Modification is allowed at the peptide N-terminus
        /// </summary>
        PeptideNTerm,

        /// <summary>
        /// Modification is allowed at the peptide C-terminus
        /// </summary>
        PeptideCTerm,

        /// <summary>
        /// Modification is allowed at the protein N-terminus
        /// </summary>
        ProteinNTerm,

        /// <summary>
        /// Modification is allowed at the protein C-terminus
        /// </summary>
        ProteinCTerm,
    }
}
