using System;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Helper class for determining what resources the system has that we can use, and for trying to not overload the system.
    /// </summary>
    public static class ParallelizationUtils
    {
        /// <summary>
        /// Number of physical cores in the system
        /// </summary>
        public static int NumPhysicalCores { get; private set; }

        /// <summary>
        /// Number of physical processors (sockets used) in the system
        /// </summary>
        public static int NumPhysicalProcessors { get; private set; }

        /// <summary>
        /// The number of logical cores in the system (includes hyperthreading cores)
        /// </summary>
        public static int NumLogicalCores { get; private set; }

        static ParallelizationUtils()
        {
            GetCoreData();
        }

        private static void GetCoreData()
        {
            NumPhysicalCores = PRISM.SystemInfo.GetCoreCount();
            NumPhysicalProcessors = PRISM.SystemInfo.GetProcessorPackageCount();
            NumLogicalCores = PRISM.SystemInfo.GetLogicalCoreCount();
        }
    }
}
