namespace InformedProteomics.Backend.Results
{
    using PSI_Interface.IdentData;

    public class ResultFileReader
    {
        /// <summary>
        /// File extensions supported by the results reader.
        /// </summary>
        public static readonly string[] SupportedResultsFiles = { ".mzid", ".mzid.gz", "_IcTda.tsv" };

        public static SimpleMZIdentMLReader.SimpleMZIdentMLData ReadResultFile(string filePath)
        {
            SimpleMZIdentMLReader.SimpleMZIdentMLData results = null;

            if (filePath.EndsWith(".mzid") || filePath.EndsWith(".mzid.gz"))
            {
                var mzidReader = new SimpleMZIdentMLReader();
                results = mzidReader.Read(filePath);
            }
            else if (filePath.EndsWith("_IcTda.tsv"))
            {
                results = DatabaseSearchResultData.ReadResultsFromFileToMzIdData(filePath);
            }

            return results;
        }
    }
}
