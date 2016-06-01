using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Database;
using PNNLOmics.Utilities;

namespace InformedProteomics.Backend.Utils
{
    public class FdrCalculator
    {
        private IList<string> _headers;
        private string[] _results;
        private readonly bool _multiplePeptidesPerScan;

        private bool _hasEvalueColumn;

        public int NumPsms { get; private set; }
        public int NumPeptides { get; private set; }

        public string ErrorMessage { get; private set; }

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

            // Add "Qvalue"
            if (!CalculateQValues(targetResultFilePath, decoyResultFilePath))
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculateQValues returned false";
                return;
            }

            // Add "PepQvalue"

            if (!CalculatePepQValues(targetResultFilePath, decoyResultFilePath))
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "CalculatePepQValues returned false";
            }
        }

        public bool HasError()
        {
            return !string.IsNullOrWhiteSpace(ErrorMessage);
        }

        public void WriteTo(string outputFilePath, bool includeDecoy = false)
        {
            var sequenceIndex = _headers.IndexOf("ProteinName");

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine(string.Join("\t", _headers));
                foreach (var r in _results)
                {
                    if (!includeDecoy && r.Split('\t')[sequenceIndex].StartsWith(FastaDatabase.DecoyProteinPrefix))
                    {
                        continue;
                    }
                    writer.WriteLine(r);
                }
            }
        }

        private bool CalculateQValues(string targetResultFilePath, string decoyResultFilePath)
        {
            string[] concatenated;
            int scoreIndex;
            int rawScoreIndex;

            var success = ReadTargetAndDecoy("QValues", targetResultFilePath, decoyResultFilePath, out concatenated, out scoreIndex, out rawScoreIndex);

            if (!success)
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "ReadTargetAndDecoy returned false in CalculateQValues";

                return false;
            }

            int scanNumIndex;
            int proteinIndex;

            if (!GetColumnIndex("QValues", "Scan", out scanNumIndex))
                return false;

            if (!GetColumnIndex("QValues", "ProteinName", out proteinIndex))
                return false;

            var distinctSorted = (_hasEvalueColumn)
                ? concatenated.OrderBy(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                    .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                    .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                    .Select(grp => grp.First())
                    .ToArray()
                : concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                    .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                    .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                    .Select(grp => grp.First())
                    .ToArray();

            NumPsms = 0;

            // Calculate q values
            _headers = _headers.Concat(new[] {"QValue"}).ToArray();
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                var columns = row.Split('\t');
                var protein = columns[proteinIndex];
                if (protein.StartsWith(FastaDatabase.DecoyProteinPrefix))
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

            _results = distinctSorted.Select((r, i) => r + "\t" + StringUtilities.DblToString(qValue[i], 7)).ToArray();

            return true;
        }

        private bool CalculatePepQValues(string targetResultFilePath, string decoyResultFilePath)
        {
            string[] concatenated;
            int scoreIndex;
            int rawScoreIndex;

            var success = ReadTargetAndDecoy("PepQValues", targetResultFilePath, decoyResultFilePath, out concatenated, out scoreIndex, out rawScoreIndex);

            if (!success)
            {
                if (string.IsNullOrWhiteSpace(ErrorMessage))
                    ErrorMessage = "ReadTargetAndDecoy returned false in CalculateQValues";

                return false;
            }

            int sequenceIndex;
            int scanNumIndex;
            int preIndex;
            int postIndex;
            int proteinIndex;

            if (!GetColumnIndex("PepQValues", "Sequence", out sequenceIndex))
                return false;

            if (!GetColumnIndex("PepQValues", "Scan", out scanNumIndex))
                return false;

            if (!GetColumnIndex("PepQValues", "Pre", out preIndex))
                return false;

            if (!GetColumnIndex("PepQValues", "Post", out postIndex))
                return false;

            if (!GetColumnIndex("PepQValues", "ProteinName", out proteinIndex))
                return false;

            string[] distinctSorted;

            if (_hasEvalueColumn)
            {
                distinctSorted = !_multiplePeptidesPerScan ?
                                concatenated.OrderBy(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                                .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                                    .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                                    .Select(grp => grp.First())
                                    .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                                    .Select(grp => grp.First())
                                    .ToArray() :
                                concatenated.OrderBy(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                                    .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                                    .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                                    .Select(grp => grp.First())
                                    .ToArray();
            }
            else
            {
                distinctSorted = !_multiplePeptidesPerScan ?
                                               concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                                               .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                                                   .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                                                   .Select(grp => grp.First())
                                                   .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                                                   .Select(grp => grp.First())
                                                   .ToArray() :
                                               concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                                                   .ThenByDescending(r => Convert.ToDouble(r.Split('\t')[rawScoreIndex]))
                                                   .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                                                   .Select(grp => grp.First())
                                                   .ToArray();
            }


            // Calculate q values
            _headers = _headers.Concat(new[] { "PepQValue" }).ToArray();
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            var peptide = new string[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                var columns = row.Split('\t');
                var protein = columns[proteinIndex];
                var annotation = columns[preIndex] + "." + columns[sequenceIndex] + "." + columns[postIndex];
                if (protein.StartsWith(FastaDatabase.DecoyProteinPrefix))
                    numDecoy++;
                else
                    numTarget++;

                fdr[i] = Math.Min(numDecoy / (double)numTarget, 1.0);
                peptide[i] = annotation;
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

            for (var i = 0; i < _results.Length; i++)
            {
                var columns = _results[i].Split('\t');
                var annotation = columns[preIndex] + "." + columns[sequenceIndex] + "." + columns[postIndex];
                _results[i] = _results[i] + "\t" + StringUtilities.DblToString(annotationToPepQValue[annotation], 7);
            }

            return true;
        }


        private bool GetColumnIndex(string targetStatistic, string columnName, out int columnIndex)
        {
            columnIndex = _headers.IndexOf(columnName);

            if (columnIndex >= 0)
            {
                return true;
            }

            ErrorMessage = "Cannot compute " + targetStatistic + "; " + columnName + " column is missing";
            return false;
        }

        private bool ReadTargetAndDecoy(
            string targetStatistic,
            string targetResultFilePath,
            string decoyResultFilePath,
            out string[] concatenated,
            out int scoreIndex,
            out int rawScoreIndex)
        {
            scoreIndex = -1;
            rawScoreIndex = -1;
            concatenated = new string[0];

            var errorBase = "Cannot compute " + targetStatistic + "; ";

            var targetData = File.ReadAllLines(targetResultFilePath);
            var decoyData = File.ReadAllLines(decoyResultFilePath);

            if (targetData.Length < 1)
            {
                ErrorMessage = errorBase + "target file is empty";
                return false;
            }

            if (decoyData.Length < 1)
            {
                ErrorMessage = errorBase + "decoy file is empty";
                return false;
            }

            var targetHeaders = targetData[0].Split('\t');
            var decoyHeaders = decoyData[0].Split('\t');

            if (targetHeaders.Length != decoyHeaders.Length)
            {
                ErrorMessage = errorBase + "header count doesn't match between target and decoy file";
                return false;
            }

            if (!targetHeaders.SequenceEqual(decoyHeaders))
            {
                ErrorMessage = errorBase + "Headers don't match between target and decoy file";
                return false;
            }

            if (_headers == null)
                _headers = targetHeaders;

            concatenated = decoyData.Skip(1).Concat(targetData.Skip(1)).ToArray();
            if (concatenated.Length == 0)
            {
                // NOTE: The DMS Analysis Manager looks for the text "No results found"
                // thus, do not change this message
                ErrorMessage = "No results found; cannot compute " + targetStatistic;
                return false;
            }

            scoreIndex = _headers.IndexOf("EValue");
            if (scoreIndex < 0)
            {
                scoreIndex = _headers.IndexOf("#MatchedFragments");
                _hasEvalueColumn = false;
            }
            else
            {
                _hasEvalueColumn = true;
            }

            if (scoreIndex < 0)
            {
                ErrorMessage = errorBase + "EValue(Score) column is missing";
                return false;
            }

            rawScoreIndex = _headers.IndexOf("Probability");
            if (rawScoreIndex < 0)
            {
                rawScoreIndex = _headers.IndexOf("#MatchedFragments");
            }

            if (rawScoreIndex < 0)
            {
                ErrorMessage = errorBase + "#MatchedFragments column is missing";
                return false;
            }

            return true;
        }

    }
}
