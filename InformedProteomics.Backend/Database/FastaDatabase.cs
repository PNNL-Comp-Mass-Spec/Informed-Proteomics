using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using ProteinFileReader;

namespace InformedProteomics.Backend.Database
{
    public class FastaDatabase
    {
        public const string SeqFileExtension = ".icseq";
        public const string AnnotationFileExtension = ".icanno";
        public const string DecoyDatabaseFileExtension = ".icdecoy.fasta";
        public const string DecoyProteinPrefix = "XXX";
        public const byte Delimiter = (byte)'_';
        public const byte LastCharacter = (byte)'~';

        private static readonly ASCIIEncoding Encoding = new ASCIIEncoding();

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

            if (!File.Exists(_seqFilePath) || !File.Exists(_annoFilePath) || !CheckHashCodeBinaryFile(_seqFilePath, _lastWriteTimeHash) ||
                !CheckHashCodeTextFile(_annoFilePath, _lastWriteTimeHash))
            {
                Console.Write("Generating " + _seqFilePath + " and " + _annoFilePath + "...");
                GenerateMetaFiles();
                Console.WriteLine("\tDone.");
            }

            IsDecoy = isDecoy;
        }

        public bool IsDecoy { get; private set; }

        public FastaDatabase Decoy(Enzyme enzyme)
        {
            if (IsDecoy)
            {
                throw new InvalidOperationException("Already a decoy database");
            }

            var decoyDatabasePath = GetDecoyDatabasePath(enzyme);
            if (!File.Exists(decoyDatabasePath)) CreateDecoyDatabase(enzyme);
            return new FastaDatabase(decoyDatabasePath, true);
        }

        public void CreateDecoyDatabase(Enzyme enzyme)
        {
            var reader = new FastaFileReader();
            if (!reader.OpenFile(_databaseFilePath))
                return;

            var decoyDatabaseFileName = GetDecoyDatabasePath(enzyme);

            Console.WriteLine("Creating " + decoyDatabaseFileName);
            using (var decoyWriter = new StreamWriter(decoyDatabaseFileName))
            {
                while (reader.ReadNextProteinEntry())
                {
                    var name = reader.ProteinName;
                    var description = reader.ProteinDescription;
                    var sequence = reader.ProteinSequence;

                    decoyWriter.WriteLine(">{0}_{1} {2}", DecoyProteinPrefix, name, description);

                    var decoySequence = new StringBuilder();
                    
                    for (var i = sequence.Length - 1; i >= 0; i--)
                    {
                        var residue = sequence[i];
                        if (enzyme.Residues.Length > 0 && enzyme.IsCleavable(residue) && decoySequence.Length > 0)
                        {
                            var residueToBeReplaced = decoySequence[decoySequence.Length-1];
                            decoySequence.Remove(decoySequence.Length-1, 1);
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
            Console.WriteLine(_sequence == null ? "Sequence is null!" : Encoding.GetString(_sequence));
        }

        public string GetDecoyDatabasePath(Enzyme enzyme)
        {
            string newExtenstion;
            if (enzyme.Residues.Length > 0) newExtenstion = ".icdecoy." + new string(enzyme.Residues) + ".fasta";
            else newExtenstion = ".icdecoy.fasta";

            return Path.ChangeExtension(_databaseFilePath, newExtenstion);
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
                fs.Seek(-sizeof(int), SeekOrigin.End);
                var lastWriteTimeHash = reader.ReadInt32();
                if (lastWriteTimeHash == code)
                    return true;
            }
            return false;
        }

        static internal bool CheckHashCodeTextFile(string filePath, int code)
        {
            var lastLine = File.ReadLines(filePath).Last(); // TODO: this is not efficient for big files
            var lastWriteTimeHash = Convert.ToInt32(lastLine);
            if (lastWriteTimeHash == code) return true;

            return false;
        }

        private readonly string _databaseFilePath;
        private readonly string _seqFilePath;
        private readonly string _annoFilePath;
        private readonly int _lastWriteTimeHash;

        private IDictionary<long, string> _names;
        private IDictionary<long, string> _descriptions;

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
                    seqWriter.Write(Encoding.GetBytes(sequence));
                    annoWriter.WriteLine(offset + ":" + name + ":" + description);
                    offset += sequence.Length;
                }

                var hashCode = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

                seqWriter.Write((char)Delimiter);
                seqWriter.Write(hashCode);
                annoWriter.Write(hashCode);

                reader.CloseFile();
            }
        }

        private bool ReadSeqFile()
        {
            using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
            {
                _sequence = new byte[fileStream.Length - sizeof(int) + 1];
                fileStream.Read(_sequence, 0, _sequence.Length);
                _sequence[_sequence.Length - 1] = LastCharacter;
            }

            return true;
        }

        //public long NumAminoAcids()
        //{
        //    using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
        //    {
        //    }
        //}

        private bool ReadAnnnoFile()
        {
            _names = new Dictionary<long, string>();
            _descriptions = new Dictionary<long, string>();

            using (var reader = new StreamReader(_annoFilePath))
            {
                string s;
                while((s=reader.ReadLine()) != null)
                {
                    var token = s.Split(':');
                    if (token.Length != 3)
                        break;
                    var offset = long.Parse(token[0]);
                    _names.Add(offset,token[1]);
                    _descriptions.Add(offset, token[2]);
                }
            }

            return true;
        }

    }
}
