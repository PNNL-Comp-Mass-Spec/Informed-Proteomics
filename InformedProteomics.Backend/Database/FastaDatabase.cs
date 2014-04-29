using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;
using ProteinFileReader;

namespace InformedProteomics.Backend.Database
{
    public class FastaDatabase
    {
        public const int FileFormatId = 175;
        public const string SeqFileExtension = ".icseq";
        public const string AnnotationFileExtension = ".icanno";
        public const string DecoyDatabaseFileExtension = ".icdecoy.fasta";
        public const string DecoyProteinPrefix = "XXX";
        public const byte Delimiter = (byte)'_';
        public const byte LastCharacter = (byte)'~';
        public const char AnnotationDelimiter = '/';

        // For suffled decoys
        public const int NumMutations = 3;

        public static readonly ASCIIEncoding Encoding = new ASCIIEncoding();

        public FastaDatabase(string databaseFilePath, bool isDecoy = false)
        {
            if (!File.Exists(databaseFilePath))
                throw new FileNotFoundException("File not found: " + databaseFilePath);
            if (string.Compare(Path.GetExtension(databaseFilePath), ".fasta", StringComparison.OrdinalIgnoreCase) != 0)
                throw new FormatException("Not a fasta file: " + databaseFilePath);

            _databaseFilePath = databaseFilePath;
            _lastWriteTimeHash = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

            var databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);
            _seqFilePath = databaseFilePathNoExt + SeqFileExtension;
            _annoFilePath = databaseFilePathNoExt + AnnotationFileExtension;

            if (!File.Exists(_seqFilePath) 
                || !File.Exists(_annoFilePath) 
                || !CheckHashCodeBinaryFile(_seqFilePath, _lastWriteTimeHash) 
                || !CheckHashCodeTextFile(_annoFilePath, _lastWriteTimeHash))
            {
                Console.Write("Generating " + _seqFilePath + " and " + _annoFilePath + "...");
                GenerateMetaFiles();
                Console.WriteLine("\tDone.");
            }

            IsDecoy = isDecoy;
        }

        public bool IsDecoy { get; private set; }

        public FastaDatabase Decoy(Enzyme enzyme, bool shuffle = false)
        {
            if (IsDecoy)
            {
                throw new InvalidOperationException("Already a decoy database");
            }

            var decoyDatabasePath = GetDecoyDatabasePath(enzyme, shuffle);
            if (!File.Exists(decoyDatabasePath)) CreateDecoyDatabase(enzyme, shuffle);
            return new FastaDatabase(decoyDatabasePath, true);
        }

        public void CreateDecoyDatabase(Enzyme enzyme, bool shuffle)
        {
            var reader = new FastaFileReader();
            if (!reader.OpenFile(_databaseFilePath))
                return;

            var decoyDatabaseFileName = GetDecoyDatabasePath(enzyme, shuffle);

            Console.WriteLine("Creating " + decoyDatabaseFileName);
            using (var decoyWriter = new StreamWriter(decoyDatabaseFileName))
            {
                while (reader.ReadNextProteinEntry())
                {
                    var name = reader.ProteinName;
                    var description = reader.ProteinDescription;
                    var sequence = reader.ProteinSequence;

                    decoyWriter.WriteLine(">{0}_{1} {2}", DecoyProteinPrefix, name, description);

                    if (!shuffle)   // Reverse
                    {
                        var decoySequence = new StringBuilder();
                        for (var i = sequence.Length - 1; i >= 0; i--)
                        {
                            var residue = sequence[i];
                            if (enzyme != null && enzyme.Residues.Length > 0 && enzyme.IsCleavable(residue) && decoySequence.Length > 0)
                            {
                                var residueToBeReplaced = decoySequence[decoySequence.Length - 1];
                                decoySequence.Remove(decoySequence.Length - 1, 1);
                                decoySequence.Append(residue);
                                decoySequence.Append(residueToBeReplaced);
                            }
                            else
                            {
                                decoySequence.Append(residue);
                            }
                        }
                        decoyWriter.WriteLine(decoySequence);
                    }
                    else
                    {
                        decoyWriter.WriteLine(SimpleStringProcessing.Mutate(SimpleStringProcessing.Shuffle(sequence), NumMutations));
                    }
                }
                reader.CloseFile();
            }
        }

        public IEnumerable<byte> Characters()
        {
            if (_sequence != null)
            {
                for (var i = 0; i < _sequence.Length - 1; i++ )
                {
                    yield return _sequence[i];
                }
            }
            else
            {
                const int bufferSize = 1 << 16;
                var buffer = new byte[bufferSize];
                var count = bufferSize;
                var numBytesRead = 0;

                using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
                {
                    var numBytesToRead = fileStream.Length - sizeof (int);
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
        }

        //public IEnumerable<byte> Characters(int index, int numCharacters)
        //{
        //    if (_sequence == null) throw new SystemException("Fasta sequence must be loaded!");
        //    for (var i = index; i < index+numCharacters && i < _sequence.Length - 1; i++)
        //    {
        //        yield return _sequence[i];
        //    }
        //}

        public void Read()
        {
            if (!ReadSeqFile())
                throw new FormatException("Error while reading " + _seqFilePath);
            if (!ReadAnnnoFile())
                throw new FormatException("Error while reading " + _annoFilePath);
        }
            
        public string GetFastaFilePath()
        {
            return _databaseFilePath;
        }

        public byte[] GetSequence()
        {
            if(_sequence == null)   Read();
            return _sequence;
        }

        public void PrintSequence()
        {
            Console.WriteLine(_sequence == null ? "Annotation is null!" : Encoding.GetString(_sequence));
        }

        public string GetDecoyDatabasePath(Enzyme enzyme, bool shuffle = false)
        {
            string newExtenstion;

            if (!shuffle)   // reverse
            {
                if (enzyme != null && enzyme.Residues.Length > 0) newExtenstion = ".icdecoy." + new string(enzyme.Residues) + ".fasta";
                else newExtenstion = ".icdecoy.fasta";
            }
            else
            {
                newExtenstion = ".icsfldecoy.fasta";
            }

            return Path.ChangeExtension(_databaseFilePath, newExtenstion);
        }

        public string GetProteinName(long offset)
        {
            var offsetKey = GetOffsetKey(offset);
            return _names[offsetKey];
        }

        public string GetProteinDescription(long offset)
        {
            var offsetKey = GetOffsetKey(offset);
            return _descriptions[offsetKey];
        }

        public int GetProteinLength(string name)
        {
            int length;
            if (_nameToLength.TryGetValue(name, out length)) return length;
            return -1;
        }

        public string GetProteinSequence(string name)
        {
            long offset;
            if (!_nameToOffset.TryGetValue(name, out offset)) return null;

            var length = _nameToLength[name];
            return Encoding.GetString(_sequence, (int)(offset + 1), length);
        }

        public int GetZeroBasedPositionInProtein(long offset)
        {
            return (int)(offset - GetOffsetKey(offset));
        }

        internal int GetLastWriteTimeHash()
        {
            return _lastWriteTimeHash;
        }

        static internal bool CheckHashCodeBinaryFile(string filePath, int code)
        {
            var fs = File.OpenRead(filePath);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-2*sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId != FileFormatId) return false;

                var lastWriteTimeHash = reader.ReadInt32();
                if (lastWriteTimeHash == code) return true;
            }
            return false;
        }

        static internal bool CheckHashCodeTextFile(string filePath, int code)
        {
            var lastLine = File.ReadLines(filePath).Last(); // TODO: this is not efficient for big files

            var token = lastLine.Split(AnnotationDelimiter);
            if (token.Length != 2) return false;

            var fileFormatId = Convert.ToInt32(token[0]);
            if (fileFormatId != FileFormatId) return false;

            var lastWriteTimeHash = Convert.ToInt32(token[1]);
            if (lastWriteTimeHash == code) return true;

            return false;
        }

        private readonly string _databaseFilePath;
        private readonly string _seqFilePath;
        private readonly string _annoFilePath;
        private readonly int _lastWriteTimeHash;

        private List<long> _offsetList;
        private IDictionary<string, int> _nameToLength; // name -> length
        private IDictionary<long, string> _names;   // offsetKey -> name
        private IDictionary<long, string> _descriptions;    // offsetKey -> description
        private IDictionary<string, long> _nameToOffset;    // name -> offsetKey

        private byte[] _sequence;

        private void GenerateMetaFiles()
        {
            if(File.Exists(_seqFilePath))
                File.Delete(_seqFilePath);

            if(File.Exists(_annoFilePath))
                File.Delete(_annoFilePath);

            using (var seqWriter = new BinaryWriter(File.Open(_seqFilePath, FileMode.CreateNew)))
            using (var annoWriter = new StreamWriter(_annoFilePath))
            {
                // Read
                var reader = new FastaFileReader();
                if (!reader.OpenFile(_databaseFilePath))
                    return;
                long offset = 0;
                while (reader.ReadNextProteinEntry())
                {
                    var name = reader.ProteinName;
                    var description = reader.ProteinDescription;
                    var sequence = (char)Delimiter + reader.ProteinSequence;
                    var length = reader.ProteinSequence.Length;
                    seqWriter.Write(Encoding.GetBytes(sequence));
                    annoWriter.WriteLine("{0}{1}{2}{3}{4}{5}{6}", 
                        offset,
                        AnnotationDelimiter,
                        length,
                        AnnotationDelimiter,
                        name,
                        AnnotationDelimiter,
                        description);
                    offset += sequence.Length;
                }

                seqWriter.Write((char)Delimiter);

                // write file format Id
                seqWriter.Write(FileFormatId);

                // writer last write time hash
                var hashCode = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

                seqWriter.Write(hashCode);
                annoWriter.Write("{0}{1}{2}", FileFormatId, AnnotationDelimiter, hashCode);

                reader.CloseFile();
            }
        }

        private bool ReadSeqFile()
        {
            using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
            {
                _sequence = new byte[fileStream.Length - 2*sizeof(int) + 1];
                fileStream.Read(_sequence, 0, _sequence.Length);
                _sequence[_sequence.Length - 1] = LastCharacter;
            }

            return true;
        }

        private bool ReadAnnnoFile()
        {
            _offsetList = new List<long>();
            _nameToLength = new Dictionary<string, int>();
            _names = new Dictionary<long, string>();
            _descriptions = new Dictionary<long, string>();
            _nameToOffset = new Dictionary<string, long>();

            using (var reader = new StreamReader(_annoFilePath))
            {
                string s;
                while((s=reader.ReadLine()) != null)
                {
                    var token = s.Split(AnnotationDelimiter);
                    if (token.Length < 4)
                        break;
                    var offset = long.Parse(token[0]);
                    _offsetList.Add(offset);
                    var length = int.Parse(token[1]);
                    var name = token[2];
                    if(_nameToLength.ContainsKey(name)) Console.WriteLine("Duplicate Name: {0}", name);
                    _nameToLength.Add(name, length);
                    _names.Add(offset,name);
                    _nameToOffset.Add(name, offset);
                    _descriptions.Add(offset, token[3]);
                }
            }

            return true;
        }

        private long GetOffsetKey(long offset)
        {
            var index = _offsetList.BinarySearch(offset);
            return index >= 0 ? _offsetList[index] : _offsetList[~index - 1];
        }
    }
}
