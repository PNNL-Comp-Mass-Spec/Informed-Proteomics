using System;

namespace InformedProteomics.Backend.Database
{
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
        Both = 3,
    }
}