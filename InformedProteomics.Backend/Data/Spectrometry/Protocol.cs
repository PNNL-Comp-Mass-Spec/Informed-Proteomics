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
        /// ITraq Label reagents
        /// </summary>
        Itraq,

        /// <summary>
        /// ITraq Label reagents with phosphorylation
        /// </summary>
        ItraqPhospho,

        /// <summary>
        /// TMT Label reagents
        /// </summary>
        Tmt
    }
}
