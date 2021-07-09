using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.SearchResults
{
    /// <summary>
    /// Extension methods for the namespace
    /// </summary>
    public static class Extensions
    {
        // Ignore Spelling: namespace

        /// <summary>
        /// Get peptides that pass the specified threshold
        /// </summary>
        /// <param name="results"></param>
        /// <param name="pepQValueThreshold"></param>
        /// <returns>List of peptides</returns>
        public static ISet<string> GetPeptides(this IEnumerable<DatabaseSearchResultData> results, double pepQValueThreshold)
        {
            return results.GetPeptidesAboveThreshold(x => x.PepQValue, pepQValueThreshold);
        }

        /// <summary>
        /// Get peptides that pass the specified threshold
        /// </summary>
        /// <param name="results"></param>
        /// <param name="qValueThreshold"></param>
        /// <returns>List of filter-passing peptides</returns>
        public static ISet<string> GetPeptidesAboveQValueThreshold(this IEnumerable<DatabaseSearchResultData> results, double qValueThreshold)
        {
            return results.GetPeptidesAboveThreshold(x => x.QValue, qValueThreshold);
        }

        /// <summary>
        /// Gets the peptides where the specified field passes the provided threshold
        /// </summary>
        /// <param name="results"></param>
        /// <param name="fieldSelector">A function, e.g. x => x.QValue</param>
        /// <param name="threshold"></param>
        /// <returns>List of filter-passing peptides</returns>
        private static ISet<string> GetPeptidesAboveThreshold(this IEnumerable<DatabaseSearchResultData> results, Func<DatabaseSearchResultData, double> fieldSelector, double threshold)
        {
            var peptideSet = new HashSet<string>();

            foreach (var result in results)
            {
                if (fieldSelector(result) <= threshold && !result.ProteinName.StartsWith("XXX"))
                {
                    peptideSet.Add(result.Sequence);
                }
            }

            return peptideSet;
        }
    }
}
