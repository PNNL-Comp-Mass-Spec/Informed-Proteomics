using System;

namespace InformedProteomics.Backend.Utils
{
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
        /// If the physical core count was retrieved via WMI, and is thereby guaranteed to be correct
        /// </summary>
        public static bool PhysCoreCountFromWmi { get; private set; }

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
            // Try to get the number of physical cores in the system - requires System.Management.dll and a WMI query, but the performance penalty for
            // using the number of logical processors in a hyperthreaded system is significant, and worse than the penalty for using fewer than all physical cores.
            NumPhysicalCores = 0;
            NumPhysicalProcessors = 0;
            PhysCoreCountFromWmi = false;

            NumLogicalCores = System.Environment.ProcessorCount;

            try
            {
                foreach (var item in new System.Management.ManagementObjectSearcher("Select NumberOfCores from Win32_Processor").Get())
                {
                    NumPhysicalProcessors++;
                    NumPhysicalCores += int.Parse(item["NumberOfCores"].ToString());
                }
                PhysCoreCountFromWmi = true;
                //Console.WriteLine("Number Of Cores: {0}", coreCount);
            }
            catch (Exception)
            {
                // Use the logical processor count, divided by 2 to avoid the greater performance penalty of over-threading.
                NumPhysicalCores = (int)(Math.Ceiling(System.Environment.ProcessorCount / 2.0));
                PhysCoreCountFromWmi = false;
            }
        }
    }
}
