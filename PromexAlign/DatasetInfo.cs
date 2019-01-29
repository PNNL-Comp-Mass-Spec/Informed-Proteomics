using System.Collections.Generic;

using System.IO;

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
        /// Parse a dataset info file to get the datasets to run on.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static List<DatasetInfo> ParseDatasetInfoFile(string filePath)
        {
            var datasets = new List<DatasetInfo>();
            var datasetNumber = 0;

            foreach (var line in File.ReadLines(filePath))
            {
                var parts = line.Split('\t');
                if (parts.Length < 3)
                {
                    continue;
                }

                datasetNumber++;

                var msPathFinderIdFilePath = string.Empty;
                if (parts.Length > 3)
                {
                    msPathFinderIdFilePath = parts[3];
                }

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
