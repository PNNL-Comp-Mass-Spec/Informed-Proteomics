using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    public enum DataFileFormat
    {
        Mgf,
        IcBottomUp,
        Dia
    };
    public static class DataFileReaderFactory
    {
        public static IDataFileReader Create(DataFileFormat format, string annotations, bool decoy, LazyLcMsRun lcms = null)
        {
            IDataFileReader reader;
            switch (format)
            {
                case DataFileFormat.Mgf:
                    reader = new MgfReader(annotations, decoy);
                    break;
                case DataFileFormat.IcBottomUp:
                    reader = new IcBottomUpTsvReader(annotations, lcms, decoy);
                    break;
                case DataFileFormat.Dia:
                    reader = new DiaTsvReader(annotations, lcms, decoy);
                    break;
                default:
                    reader = null;
                    break;
            }
            return reader;
        }
    }
}
