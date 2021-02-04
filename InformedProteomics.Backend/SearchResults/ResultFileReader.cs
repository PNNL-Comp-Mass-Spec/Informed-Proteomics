using PSI_Interface.IdentData;

namespace InformedProteomics.Backend.SearchResults
{
    /// <summary>
    /// Reader factory for result files
    /// </summary>
    public static class ResultFileReader
    {
        /// <summary>
        /// File extensions supported by the results reader.
        /// These extensions need to be lowercase
        /// </summary>
        public static readonly string[] SupportedResultsFiles = { ".mzid", ".mzid.gz", "_ictda.tsv" };

        /// <summary>
        /// Read the file at path <paramref name="filePath"/>.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static SimpleMZIdentMLReader.SimpleMZIdentMLData ReadResultFile(string filePath)
        {
            SimpleMZIdentMLReader.SimpleMZIdentMLData results = null;

            var lowerFilePath = filePath.ToLower();

            if (lowerFilePath.EndsWith(".mzid") || lowerFilePath.EndsWith(".mzid.gz"))
            {
                var mzidReader = new SimpleMZIdentMLReader();
                results = mzidReader.Read(filePath);
            }
            else if (lowerFilePath.EndsWith("_ictda.tsv"))
            {
                results = DatabaseSearchResultData.ReadResultsFromFileToMzIdData(filePath);
            }

            return results;
        }
    }
}
