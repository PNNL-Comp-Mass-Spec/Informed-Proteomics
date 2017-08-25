namespace InformedProteomics.Backend.Database
{
    public static class FastaDatabaseConstants
    {
        /// <summary>
        /// Prefix to flag decoy proteins
        /// </summary>
        public const string DecoyProteinPrefix = "XXX";

        /// <summary>
        /// Sequence delimiter in the backing files
        /// </summary>
        public const byte Delimiter = (byte)'_';

        /// <summary>
        /// Last character marker in the backing files
        /// </summary>
        public const byte LastCharacter = (byte)'~';

        /// <summary>
        /// Annotation delimiter in the backing files
        /// </summary>
        public const char AnnotationDelimiter = '/';
    }
}
