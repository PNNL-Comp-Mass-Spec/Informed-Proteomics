using System.IO;

namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Constants used in multiple places where Fasta Databases are used
    /// </summary>
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

        /// <summary>
        /// Examine the extension on the database file
        /// Return true if .fasta or .fa or .faa
        /// Otherwise, return false
        /// </summary>
        /// <param name="databaseFilePath"></param>
        /// <returns></returns>
        public static bool ValidFASTAExtension(string databaseFilePath)
        {
            if (string.IsNullOrWhiteSpace(databaseFilePath))
                return false;

            var dbExtension = Path.GetExtension(databaseFilePath).ToLower();

            return dbExtension.Equals(".fasta") || dbExtension.Equals(".fa") || dbExtension.Equals(".faa");
        }
    }
}
