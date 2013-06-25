using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ProteinFileReader;

namespace InformedProteomics.Backend.Database
{
    public class FastaDatabase
    {
        public static readonly string SeqFileExtension = ".icseq";
        public static readonly string AnnotationFileExtension = ".icanno";
        public static readonly char Delimiter = '_';
        public static readonly byte LastCharacter = (byte)'~';

        private static readonly ASCIIEncoding Encoding = new ASCIIEncoding();

        public FastaDatabase(string databaseFilePath)
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
                Console.Write("Generating " + _seqFilePath + " and " + _annoFilePath);
                GenerateMetaFiles();
                Console.WriteLine("\tDone.");
            }
        }

        public IEnumerable<char> Characters()
        {
            if (_sequence != null)
            {
                foreach (var code in _sequence)
                {
                    yield return Convert.ToChar(code);
                }
            }
            else
            {
                using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
                {
                    for (var i = 0; i < fileStream.Length - sizeof (int); i++)
                    {
                        yield return Convert.ToChar(fileStream.ReadByte());
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

        internal int GetLastWriteTimeHash()
        {
            return _lastWriteTimeHash;
        }

        static internal bool CheckHashCodeBinaryFile(string filePath, int code)
        {
            using (var fs = File.OpenRead(filePath))
            {
                fs.Seek(-sizeof(int), SeekOrigin.End);

                using (var reader = new BinaryReader(fs))
                {
                    int lastWriteTimeHash = reader.ReadInt32();
                    if (lastWriteTimeHash == code)
                        return true;
                }
            }
            return false;
        }

        static internal bool CheckHashCodeTextFile(string filePath, int code)
        {
            var lastLine = File.ReadLines(filePath).Last(); // TODO: this is not efficient for big files
            int lastWriteTimeHash = Convert.ToInt32(lastLine);
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
                    var sequence = Delimiter + reader.ProteinSequence;
                    seqWriter.Write(Encoding.GetBytes(sequence));
                    annoWriter.WriteLine(offset + ":" + name + ":" + description);
                    offset += sequence.Length;
                }

                var hashCode = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

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
