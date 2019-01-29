namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Types of MS2 ion detectors
    /// </summary>
    public enum Ms2DetectorType : byte
    {
        /// <summary>
        /// Ion trap
        /// </summary>
        IonTrap,

        /// <summary>
        /// Thermo Orbitrap
        /// </summary>
        Orbitrap,

        /// <summary>
        /// FTICR
        /// </summary>
        Fticr,

        /// <summary>
        /// Time-of-flight
        /// </summary>
        Tof
    }
}
