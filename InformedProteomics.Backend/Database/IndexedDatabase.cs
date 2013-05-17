using System.IO;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    public class IndexedDatabase
    {
        public static readonly string PermutedLongestCommonPrefixFileExtension = ".icplcp";

        private readonly byte[] _pLcp;  // permuted longest common prefix

        private readonly string _pLcpFilePath;

        public IndexedDatabase(FastaDatabase fastaDatabase)
        {
            string databaseFilePath = fastaDatabase.GetFastaFilePath();
            var databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);

            _pLcpFilePath = databaseFilePathNoExt + PermutedLongestCommonPrefixFileExtension;
            var lastWriteTimeHash = fastaDatabase.GetLastWriteTimeHash();

            if (!File.Exists(_pLcpFilePath) || !FastaDatabase.CheckHashCodeBinaryFile(_pLcpFilePath, lastWriteTimeHash))
            {
                CreatePermutedLongestCommonPrefixFile(fastaDatabase.GetSequence());
            }
        }

        private void CreatePermutedLongestCommonPrefixFile(byte[] sequence)
        {
            if (File.Exists(_pLcpFilePath))
                File.Delete(_pLcpFilePath);

            var suffixArray = new int[sequence.Length];
            SAIS.sufsort(sequence, suffixArray, sequence.Length);

        }

        private byte[] CreatePermutedLcp(byte[] sequence, int[] suffixArray)
        {
            return null;
        }
    }
}
