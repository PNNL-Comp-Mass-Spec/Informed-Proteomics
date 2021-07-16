using System;
using System.Collections.Generic;

using System.IO;
using PRISM;

namespace PromexAlign
{
    public class DatasetInfo
    {
        /// <summary>
        /// Gets or sets the label for this dataset.
        /// </summary>
        public string Label { get; set; }

        /// <summary>
        /// Gets or sets the path to the raw file.
        /// </summary>
        public string RawFilePath { get; set; }

        /// <summary>
        /// Gets the path to the ms1ft file.
        /// </summary>
        public string Ms1FtFilePath { get; set; }

        /// <summary>
        /// Gets the path to the MSPathFinder results file (_IcTda.tsv)
        /// </summary>
        public string MsPfIdFilePath { get; set; }

        /// <summary>
        /// Parse a dataset info file to get the datasets to run on
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns>List of dataset info objects</returns>
        /// <remarks>
        /// <para>
        /// The input file should be a tab-delimited file with 3 or 4 columns. The header line (with column names) is optional.
        /// </para>
        /// <para>
        /// Expected columns:
        /// Label  RawFilePath  Ms1FtFilePath  MsPathfinderIdFilePath
        /// </para>
        /// </remarks>
        public static List<DatasetInfo> ParseDatasetInfoFile(string filePath)
        {
            var fileNameShown = false;

            var datasets = new List<DatasetInfo>();
            var datasetNumber = 0;
            var rowNumber = 0;

            using var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

            while (!reader.EndOfStream)
            {
                var dataLine = reader.ReadLine();
                rowNumber++;

                if (string.IsNullOrWhiteSpace(dataLine))
                    continue;

                var parts = dataLine.Split('\t');
                if (parts.Length < 3)
                {
                    if (!fileNameShown)
                    {
                        Console.WriteLine("Loading filenames from file " + filePath);
                        fileNameShown = true;
                    }

                    ConsoleMsgUtils.ShowWarning("Skipping row {0} since it has fewer than 3 columns", rowNumber);
                    continue;
                }

                if (parts[0].Equals("Label", StringComparison.OrdinalIgnoreCase) &&
                    parts[1].Equals("RawFilePath", StringComparison.OrdinalIgnoreCase))
                {
                    // Header line; skip it
                    continue;
                }

                datasetNumber++;

                var msPathFinderIdFilePath = parts.Length > 3 ? parts[3] : string.Empty;

                var datasetInfo = new DatasetInfo
                {
                    Label = parts[0],
                    RawFilePath = parts[1],
                    Ms1FtFilePath = parts[2],
                    MsPfIdFilePath = msPathFinderIdFilePath
                };

                if (string.IsNullOrWhiteSpace(datasetInfo.Label))
                {
                    datasetInfo.Label = "Dataset_" + datasetNumber;
                }

                datasets.Add(datasetInfo);
            }

            return datasets;
        }
    }
}
