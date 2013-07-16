using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    public class IndexedDatabase
    {
        public static readonly string PermutedLongestCommonPrefixFileExtension = ".icplcp";

        private readonly string _pLcpFilePath;

        private readonly FastaDatabase _fastaDatabase;
        private readonly int _maxPeptideLength;
        private byte[] _pLcp;        

        public IndexedDatabase(FastaDatabase fastaDatabase)
        {
            _fastaDatabase = fastaDatabase;
            var databaseFilePath = _fastaDatabase.GetFastaFilePath();
            var databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);

            _pLcpFilePath = databaseFilePathNoExt + PermutedLongestCommonPrefixFileExtension;
            var lastWriteTimeHash = _fastaDatabase.GetLastWriteTimeHash();

            if (!File.Exists(_pLcpFilePath) || !FastaDatabase.CheckHashCodeBinaryFile(_pLcpFilePath, lastWriteTimeHash))
            {
                Console.Write("Generating " + _pLcpFilePath + "...");
                CreatePermutedLongestCommonPrefixFile();
                Console.WriteLine("\tDone.");
            }

            _maxPeptideLength = 30;
        }

        public void Read()
        {
            using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
            {
                _pLcp = new byte[fileStream.Length - sizeof(int)];
                fileStream.Read(_pLcp, 0, _pLcp.Length);
            }
        }

        public IEnumerable<byte> PLcps()
        {
            if (_pLcp == null)
            {
                const int bufferSize = 1 << 16;
                var buffer = new byte[bufferSize];
                var count = bufferSize;
                var numBytesRead = 0;

                using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
                {
                    var numBytesToRead = fileStream.Length - sizeof(int);
                    while (count > 0)
                    {
                        count = fileStream.Read(buffer, 0, bufferSize);
                        for (var i = 0; i < count && numBytesRead++ < numBytesToRead; i++)
                        {
                            yield return buffer[i];
                        }
                    }
                }                
            }
            else
            {
                foreach (var lcp in _pLcp) yield return lcp;
            }
        }

        public IEnumerable<Tuple<byte[], short>> SequenceLcpPairs(int minLength, int maxLength)
        {
            var curSequence = new LinkedList<byte>();
            var lcpList = new LinkedList<byte>();
            var lcpEnum = PLcps().GetEnumerator();

            foreach (var residue in _fastaDatabase.Characters())
            {
                lcpEnum.MoveNext();
                curSequence.AddLast(residue);
                lcpList.AddLast(lcpEnum.Current);

                if (curSequence.Count < maxLength) continue;

                var seqArr = new byte[curSequence.Count];
                curSequence.CopyTo(seqArr, 0);

                yield return new Tuple<byte[], short>(seqArr, lcpList.First.Value);
                curSequence.RemoveFirst();
                lcpList.RemoveFirst();
            }

            while (curSequence.Count >= minLength)
            {
                var seqArr = new byte[curSequence.Count];
                curSequence.CopyTo(seqArr, 0);
                yield return new Tuple<byte[], short>(seqArr, lcpList.First.Value);
                curSequence.RemoveFirst();
                lcpList.RemoveFirst();
            }
        }

        public long CountSequences(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme)
        {
            long numSequences = 0;

            return numSequences;
        }

        public long NumSequences(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme)
        {
            return NumSequences(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
                               enzyme.IsNTerm);
        }

        public IEnumerable<string> SequencesAsStrings(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme)
        {
            return SequencesAsStrings(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
                               enzyme.IsNTerm);
        }

        public long NumSequences(int minLength, int maxLength, int numTolerableTermini,
                                 int numMissedCleavages, char[] enzymaticResidues,
                                 bool isNTermEnzyme)
        {
            long numSequences = 0L;

            {
                var isCleavable = new bool[128];
                if (enzymaticResidues != null)
                {
                    foreach (var residue in enzymaticResidues)
                    {
                        isCleavable[(int)residue] = true;
                        isCleavable[(int)FastaDatabase.Delimiter] = true;
                    }
                }

                const string standardAAStr = "ACDEFGHIKLMNPQRSTVWY";
                var isStandardAminoAcid = new bool[128];
                foreach (var residue in standardAAStr)
                {
                    isStandardAminoAcid[(int)residue] = true;
                }

                var encoding = System.Text.Encoding.ASCII;
                // pre, peptide sequence, next
                foreach (var seqAndLcp in SequenceLcpPairs(minLength, maxLength + 2))
                {
                    var seqArr = seqAndLcp.Item1;
                    var lcp = seqAndLcp.Item2;
                    var ntt = 0;
                    var nmc = 0;

                    if (enzymaticResidues != null)
                    {
                        if (!isNTermEnzyme) // C-term enzyme 
                        {
                            if (isCleavable[seqArr[0]]) ++ntt;
                            if (ntt < numTolerableTermini - 1) continue;

                            for (var i = 1; i < seqArr.Length - 1; i++)
                            {
                                var code = seqArr[i];
                                if (!isStandardAminoAcid[code]) break;
                                if (isCleavable[code]) ++nmc;

                                if (i >= minLength && i >= lcp)
                                {
                                    if (ntt + (isCleavable[code] || i < seqArr.Length - 1 && seqArr[i + 1] == FastaDatabase.Delimiter ? 1 : 0) >= numTolerableTermini)
                                    {
                                        numSequences++;
                                    }
                                }
                                if (nmc > numMissedCleavages) break;
                            }
                        }
                        else // N-term enzyme
                        {

                        }
                    }
                    else // no enzyme
                    {
                    }
                }
            }
            return numSequences;
        }

        public IEnumerable<string> SequencesAsStrings(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, char[] enzymaticResidues,
                                                      bool isNTermEnzyme)
        {
            var isCleavable = new bool[128];
            if (enzymaticResidues != null)
            {
                foreach (var residue in enzymaticResidues)
                {
                    isCleavable[(int)residue] = true;
                    isCleavable[(int) FastaDatabase.Delimiter] = true;
                }
            }

            const string standardAAStr = "ACDEFGHIKLMNPQRSTVWY";
            var isStandardAminoAcid = new bool[128];
            foreach (var residue in standardAAStr)
            {
                isStandardAminoAcid[(int) residue] = true;
            }

            var encoding = System.Text.Encoding.ASCII;
            // pre, peptide sequence, next
            foreach (var seqAndLcp in SequenceLcpPairs(minLength, maxLength+2))
            {
                var seqArr = seqAndLcp.Item1;
                var lcp = seqAndLcp.Item2;
                var ntt = 0;
                var nmc = 0;

                if (enzymaticResidues != null)  
                {
                    if (!isNTermEnzyme) // C-term enzyme 
                    {
                        if (isCleavable[seqArr[0]]) ++ntt;
                        if (ntt < numTolerableTermini-1) continue;

                        for (var i = 1; i < seqArr.Length-1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;

                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[code] || i < seqArr.Length-1 && seqArr[i+1] == FastaDatabase.Delimiter ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return string.Format("{0}.{1}.{2}",(char)seqArr[0], encoding.GetString(seqArr, 1, i), (char)seqArr[i+1]);
                                }
                            }
                            if (nmc > numMissedCleavages) break;
                        }
                    }
                    else // N-term enzyme
                    {
                        
                    }
                }
                else // no enzyme
                {
                }
            }
        }

        private void CreatePermutedLongestCommonPrefixFile()
        {
            if (File.Exists(_pLcpFilePath))
                File.Delete(_pLcpFilePath);

            var sequence = _fastaDatabase.GetSequence();
            //Console.WriteLine("Sequence: {0}", System.Text.Encoding.ASCII.GetString(sequence));
            var suffixArray = new int[sequence.Length-1];
            SAIS.sufsort(sequence, suffixArray, sequence.Length-1);

            var prevIndex = sequence.Length - 1;

            var pLcp = new byte[suffixArray.Length];

            foreach (var index in suffixArray)
            {
                var lcp = GetLcp(sequence, prevIndex, index);
                pLcp[index] = lcp;
                //Console.WriteLine("{0}\t{1}", System.Text.Encoding.ASCII.GetString(sequence, index, sequence.Length - index - 1), lcp);
                prevIndex = index;
            }

            using (var fs = new FileStream(_pLcpFilePath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                foreach (var lcp in pLcp)
                {
                    //Console.WriteLine("LCP: {0}", lcp);
                    fs.WriteByte(lcp);
                }
                fs.Write(BitConverter.GetBytes(_fastaDatabase.GetLastWriteTimeHash()), 0, sizeof(int));
            }
        }

        private static byte GetLcp(IList<byte> sequence, int index1, int index2)
        {
            var lcp = (byte)0;

            while (sequence[index1+lcp] == sequence[index2+lcp])
            {
                ++lcp;
                if (lcp == byte.MaxValue) break;
            }

            return lcp;
        }
    }
}
