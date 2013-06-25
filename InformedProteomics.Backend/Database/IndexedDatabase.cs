using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    public class IndexedDatabase
    {
        public static readonly string PermutedLongestCommonPrefixFileExtension = ".icplcp";

        private readonly string _pLcpFilePath;

        private readonly FastaDatabase _fastaDatabase;

        private readonly int _maxPeptideLength;
        
        private LinkedList<char> _curSequence;

        public IndexedDatabase(FastaDatabase fastaDatabase)
        {
            _fastaDatabase = fastaDatabase;
            var databaseFilePath = _fastaDatabase.GetFastaFilePath();
            var databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);

            _pLcpFilePath = databaseFilePathNoExt + PermutedLongestCommonPrefixFileExtension;
            var lastWriteTimeHash = _fastaDatabase.GetLastWriteTimeHash();

            if (!File.Exists(_pLcpFilePath) || !FastaDatabase.CheckHashCodeBinaryFile(_pLcpFilePath, lastWriteTimeHash))
            {
                Console.Write("Generating " + _pLcpFilePath);
                CreatePermutedLongestCommonPrefixFile();
                Console.WriteLine("\tDone.");
            }

            _maxPeptideLength = 30;
            _curSequence = null;
        }

        public IEnumerable<Tuple<string, short>> Sequence(int maxLength)
        {
            using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
            {
                foreach (var residue in _fastaDatabase.Characters())
                {
                }
            }
            return null;
        }

        private void CreatePermutedLongestCommonPrefixFile()
        {
            if (File.Exists(_pLcpFilePath))
                File.Delete(_pLcpFilePath);

            var sequence = _fastaDatabase.GetSequence();
            //Console.WriteLine("Sequence: {0}", Encoding.ASCII.GetString(sequence));
            var suffixArray = new int[sequence.Length-1];
            SAIS.sufsort(sequence, suffixArray, sequence.Length-1);

            var prevIndex = sequence.Length - 1;

            var pLcp = new short[suffixArray.Length];

            foreach (var index in suffixArray)
            {
                var lcp = GetLcp(sequence, prevIndex, index);
                pLcp[index] = lcp;
                //Console.WriteLine("{0}\t{1}", Encoding.ASCII.GetString(sequence, index, sequence.Length - index - 1), lcp);
                prevIndex = index;
            }

            using (var fs = new FileStream(_pLcpFilePath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                using (var pLcpWriter = new BinaryWriter(fs))
                {
                    foreach(var val in pLcp)    pLcpWriter.Write(val);
                    pLcpWriter.Write(_fastaDatabase.GetLastWriteTimeHash());
                }
            }
        }

        private static short GetLcp(IList<byte> sequence, int index1, int index2)
        {
            var lcp = (short)0;

            while (sequence[index1+lcp] == sequence[index2+lcp])
            {
                ++lcp;
            }

            return lcp;
        }
    }
}
