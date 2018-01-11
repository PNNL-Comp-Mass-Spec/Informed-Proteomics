using System;
using System.ComponentModel;

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
        [Description("Target search only")]
        Target = 1,

        /// <summary>
        /// Decoy search only
        /// </summary>
        [Description("Decoy search only (shuffled database)")]
        Decoy = 2,

        /// <summary>
        /// Target and Decoy search
        /// </summary>
        [Description("Target and Decoy search")]
        Both = Target | Decoy,
    }
}