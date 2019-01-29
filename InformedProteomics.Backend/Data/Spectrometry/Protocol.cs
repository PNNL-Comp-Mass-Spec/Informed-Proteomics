namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Training protocol options
    /// </summary>
    public enum Protocol : byte
    {
        /// <summary>
        /// Standard protocol
        /// </summary>
        Standard,

        /// <summary>
        /// Phosporylation
        /// </summary>
        Phosphorylation,

        /// <summary>
        /// iTRAQ Label reagents
        /// </summary>
        Itraq,

        /// <summary>
        /// iTRAQ Label reagents with phosphorylation
        /// </summary>
        ItraqPhospho,

        /// <summary>
        /// TMT Label reagents
        /// </summary>
        Tmt
    }
}
