using System;

namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Available database search modes
    /// </summary>
    [Flags]
    public enum DatabaseSearchMode
    {
        /// <summary>
        /// Target search only
        /// </summary>
        Target = 1,

        /// <summary>
        /// Decoy search only
        /// </summary>
        Decoy = 2,

        /// <summary>
        /// Target and Decoy search
        /// </summary>
        Both = Target | Decoy,
    }
}