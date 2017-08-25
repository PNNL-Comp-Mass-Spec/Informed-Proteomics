using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Results;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Computes the False Decoy Ratio and scores for the supplied target and decoy hits
    /// </summary>
    public class FdrCalculator
    {
        private readonly bool _multiplePeptidesPerScan;
        private List<DatabaseSearchResultData> searchResults = new List<DatabaseSearchResultData>();
        private List<DatabaseSearchResultData> filteredResults = new List<DatabaseSearchResultData>();


        /// <summary>
        /// Number of PSMs with a QValue &lt; 0.01
        /// </summary>
        public int NumPsms { get; private set; }

        /// <summary>
        /// Number of peptides with a PepQValue &lt; 0.01
        /// </summary>
        public int NumPeptides { get; private set; }

        /// <summary>
        /// Error message, if FDR calculation fails
        /// </summary>
        public string ErrorMessage { get; private set; }

        /// <summary>
        /// The full list of filtered results, with FDR scores added
        /// </summary>
        public List<DatabaseSearchResultData> FilteredResults
        {
            get { return new List<DatabaseSearchResultData>(filteredResults); }
        }

        /// <summary>
        /// Instantiate the FDR calculator
        /// </summary>
        /// <param name="targetResultFilePath"></param>
        /// <param name="decoyResultFilePath"></param>
        /// <param name="multiplePeptidesPerScan"></param>
        /// <remarks>If an error occurs, ErrorMessage will be non-null</remarks>
        public FdrCalculator(string targetResultFilePath, string decoyResultFilePath, bool multiplePeptidesPerScan = false)
        {
            ErrorMessage = string.Empty;

            _multiplePeptidesPerScan = multiplePeptidesPerScan;

            // Read the data and add it to the list
            if (!ReadTargetAndDecoy(targetResultFilePath, decoyResultFilePath))
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "ReadTargetAndDecoy returned false in FdrCalculator";

                return;
            }

            // Add "Qvalue"
            if (!CalculateQValues())
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculateQValues returned false";
                return;
            }

            // Add "PepQvalue"
            if (!CalculatePepQValues())
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculatePepQValues returned false";
            }
        }

        /// <summary>
        /// Instantiate the FDR calculator
        /// </summary>
        /// <param name="targetResults"></param>
        /// <param name="decoyResults"></param>
        /// <param name="multiplePeptidesPerScan"></param>
        public FdrCalculator(List<DatabaseSearchResultData> targetResults, List<DatabaseSearchResultData> decoyResults, bool multiplePeptidesPerScan = false)
        {
            ErrorMessage = string.Empty;

            _multiplePeptidesPerScan = multiplePeptidesPerScan;

            // Add the data to the list
            if (!AddTargetAndDecoyData(targetResults, decoyResults))
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "AddTargetAndDecoy returned false in FdrCalculator";

                return;
            }

            // Add "Qvalue"
            if (!CalculateQValues())
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculateQValues returned false";
                return;
            }

            // Add "PepQvalue"
            if (!CalculatePepQValues())
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculatePepQValues returned false";
            }
        }

        /// <summary>
        /// True if there was an error calculating the FDR scores
        /// </summary>
        /// <returns></returns>
        public bool HasError()
        {
            return !string.IsNullOrWhiteSpace(ErrorMessage);
        }

        /// <summary>
        /// Write the results with the FDR data to the specified file
        /// </summary>
        /// <param name="outputFilePath"></param>
        /// <param name="includeDecoy"></param>
        public void WriteTo(string outputFilePath, bool includeDecoy = false)
        {
            //var resultsToUse = searchResults;
            var resultsToUse = filteredResults;

            IEnumerable<DatabaseSearchResultData> resultsToWrite;
            if (includeDecoy)
            {
                resultsToWrite = resultsToUse;
            }
            else
            {
                resultsToWrite = resultsToUse.Where(x => !x.ProteinName.StartsWith(FastaDatabaseConstants.DecoyProteinPrefix));
            }
            DatabaseSearchResultData.WriteResultsToFile(outputFilePath, resultsToWrite, true);
        }

        private bool CalculateQValues()
        {
            // Order by EValue, then Probability, then take only the best scoring result for each scan number
            var distinctSorted = searchResults.OrderBy(r => r.EValue)
                .ThenByDescending(r => r.Probability)
                .GroupBy(r => r.ScanNum)
                .Select(grp => grp.First())
                .ToArray();

            NumPsms = 0;

            // Calculate q values
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var result = distinctSorted[i];
                if (result.ProteinName.StartsWith(FastaDatabaseConstants.DecoyProteinPrefix))
                    numDecoy++;
                else
                    numTarget++;

                fdr[i] = Math.Min(numDecoy / (double)numTarget, 1.0);
            }

            var qValue = new double[fdr.Length];
            qValue[fdr.Length - 1] = fdr[fdr.Length - 1];
            for (var i = fdr.Length - 2; i >= 0; i--)
            {
                qValue[i] = Math.Min(qValue[i + 1], fdr[i]);
                if (qValue[i] <= 0.01)
                    ++NumPsms;
            }

            for (var i = 0; i < distinctSorted.Length; i++)
            {
                distinctSorted[i].QValue = qValue[i];
            }

            filteredResults.AddRange(distinctSorted);

            return true;
        }

        private bool CalculatePepQValues()
        {
            IEnumerable<DatabaseSearchResultData> distinctSorting = searchResults.OrderBy(r => r.EValue).ThenByDescending(r => r.Probability);
            if (!_multiplePeptidesPerScan)
            {
                distinctSorting = distinctSorting.GroupBy(r => r.ScanNum).Select(grp => grp.First());
            }
            var distinctSorted = distinctSorting.GroupBy(r => r.SequenceWithEnds).Select(grp => grp.First()).ToArray();

            // Calculate pepq values
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            var peptide = new string[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                if (row.ProteinName.StartsWith(FastaDatabaseConstants.DecoyProteinPrefix))
                    numDecoy++;
                else
                    numTarget++;

                fdr[i] = Math.Min(numDecoy / (double)numTarget, 1.0);
                peptide[i] = row.SequenceWithEnds;
            }

            var pepQValue = new double[fdr.Length];
            pepQValue[fdr.Length - 1] = fdr[fdr.Length - 1];
            for (var i = fdr.Length - 2; i >= 0; i--)
            {
                pepQValue[i] = Math.Min(pepQValue[i + 1], fdr[i]);
                if (pepQValue[i] <= 0.01)
                    ++NumPeptides;
            }

            var annotationToPepQValue = new Dictionary<string, double>();
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                annotationToPepQValue[peptide[i]] = pepQValue[i];
            }

            foreach (var item in filteredResults)
            {
                item.PepQValue = annotationToPepQValue[item.SequenceWithEnds];
            }

            return true;
        }

        private bool ReadTargetAndDecoy(string targetResultFilePath, string decoyResultFilePath)
        {
            var errorBase = "Cannot compute FDR Scores; ";

            if (!File.Exists(targetResultFilePath))
            {
                ErrorMessage = errorBase + "target file not found, " + Path.GetFileName(targetResultFilePath);
                return false;
            }

            if (!File.Exists(decoyResultFilePath))
            {
                ErrorMessage = errorBase + "decoy file not found, " + Path.GetFileName(decoyResultFilePath);
                return false;
            }

            var targetData = DatabaseSearchResultData.ReadResultsFromFile(targetResultFilePath);
            var decoyData = DatabaseSearchResultData.ReadResultsFromFile(decoyResultFilePath);

            return AddTargetAndDecoyData(targetData, decoyData);
        }

        private bool AddTargetAndDecoyData(List<DatabaseSearchResultData> targetResults, List<DatabaseSearchResultData> decoyResults)
        {
            var errorBase = "Cannot compute FDR Scores; ";
            if (targetResults == null || targetResults.Count < 1)
            {
                ErrorMessage = errorBase + "target file is empty";
                return false;
            }

            if (decoyResults == null || decoyResults.Count < 1)
            {
                ErrorMessage = errorBase + "decoy file is empty";
                return false;
            }

            searchResults = new List<DatabaseSearchResultData>();
            searchResults.AddRange(decoyResults);
            searchResults.AddRange(targetResults);

            if (searchResults.Count == 0)
            {
                // NOTE: The DMS Analysis Manager looks for the text "No results found"
                // thus, do not change this message
                ErrorMessage = "No results found; cannot compute FDR Scores";
                return false;
            }

            return true;
        }
    }
}
