using System;
using System.Reflection;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Miscellaneous utilities
    /// </summary>
    public static class Misc
    {
        // Ignore Spelling: yyyy

        /// <summary>
        /// Parse the assembly's version to try to get a build date
        /// </summary>
        /// <returns></returns>
        public static DateTime GetBuildDateFromVersion()
        {
            try
            {
                var version = Assembly.GetEntryAssembly().GetName().Version;
                return GetBuildDateFromVersion(version);
            }
            catch
            {
                // This method was likely invoked by NUnit
                var version = new Version(1, 0);
                return GetBuildDateFromVersion(version);
            }
        }

        /// <summary>
        /// Parse the supplied version to try to get a build date
        /// </summary>
        /// <param name="version"></param>
        /// <returns></returns>
        public static DateTime GetBuildDateFromVersion(Version version)
        {
            var buildDateTime = new DateTime(2000, 1, 1).Add(
                new TimeSpan(TimeSpan.TicksPerDay * version.Build + // days since 1 January 2000
                             TimeSpan.TicksPerSecond * 2 * version.Revision)); // seconds since midnight, (multiply by 2 to get original)

            return buildDateTime;
        }

        /// <summary>
        /// Parse the assembly's version to try to get a build date
        /// </summary>
        /// <returns></returns>
        public static string GetBuildDateTextFromVersion()
        {
            try
            {
                var version = Assembly.GetEntryAssembly().GetName().Version;
                return GetBuildDateTextFromVersion(version);
            }
            catch
            {
                // This method was likely invoked by NUnit
                var version = new Version(1, 0);
                return GetBuildDateTextFromVersion(version);
            }
        }

        /// <summary>
        /// Parse the supplied version to try to get a build date
        /// </summary>
        /// <param name="version"></param>
        /// <returns></returns>
        public static string GetBuildDateTextFromVersion(Version version)
        {
            var buildDateTime = GetBuildDateFromVersion(version);

            // ReSharper disable once StringLiteralTypo
            return buildDateTime.ToString("MMMM d, yyyy");
        }
    }
}
