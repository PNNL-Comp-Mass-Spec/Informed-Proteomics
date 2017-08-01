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
        CID = 0,

        /// <summary>
        /// Electron transfer dissociation
        /// </summary>
        ETD = 1,

        /// <summary>
        /// High-energy/beam-type collision induced dissociation
        /// </summary>
        HCD = 2,

        /// <summary>
        /// Electron capture dissociation
        /// </summary>
        ECD = 3,

        /// <summary>
        /// Pulsed q dissociation
        /// </summary>
        PQD = 4,

        /// <summary>
        /// Ultraviolet photo dissociation
        /// </summary>
        UVPD = 5,

        /// <summary>
        /// Unknown activation method
        /// </summary>
        Unknown = 6,
    }
}
