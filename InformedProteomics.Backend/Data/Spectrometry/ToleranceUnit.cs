namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Units used for tolerance
    /// </summary>
    public enum ToleranceUnit
    {
        /// <summary>
        /// Parts per million
        /// </summary>
        Ppm,

        /// <summary>
        /// Daltons
        /// </summary>
        Da,

        /// <summary>
        /// m/z
        /// </summary>
        Mz,

        /// <summary>
        /// Thomsons (m/z)
        /// </summary>
        [System.Obsolete("Use Mz - use of m/z is standard, Thomsons is not.")]
        Th = Mz
    }
}
