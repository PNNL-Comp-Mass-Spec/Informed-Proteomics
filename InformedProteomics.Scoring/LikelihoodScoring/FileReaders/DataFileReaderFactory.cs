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
            return format switch
            {
                DataFileFormat.Mgf => new MgfReader(annotations, decoy),
                DataFileFormat.IcBottomUp => new IcBottomUpTsvReader(annotations, lcms, decoy),
                DataFileFormat.Dia => new DiaTsvReader(annotations, lcms, decoy),
                _ => null
            };
        }
    }
}
