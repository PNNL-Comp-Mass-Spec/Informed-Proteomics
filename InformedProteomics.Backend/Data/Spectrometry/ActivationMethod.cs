namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Type of dissociation/activation supported/available in this program
    /// </summary>
    public enum ActivationMethod : byte
    {
        /// <summary>
        /// Collision-induced dissociation
        /// </summary>
        CID,

        /// <summary>
        /// Electron transfer dissociation
        /// </summary>
        ETD,

        /// <summary>
        /// High-energy/beam-type collision induced dissociation
        /// </summary>
        HCD,

        /// <summary>
        /// Electron capture dissociation
        /// </summary>
        ECD,

        /// <summary>
        /// Pulsed q dissociation
        /// </summary>
        PQD,

        /// <summary>
        /// Unknown activation method
        /// </summary>
        Unknown,
    }
}
