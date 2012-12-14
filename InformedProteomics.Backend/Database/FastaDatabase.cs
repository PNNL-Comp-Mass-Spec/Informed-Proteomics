using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using ProteinFileReader;

namespace InformedProteomics.Backend.Database
{
    public class FastaDatabase
    {
        public static readonly string SeqFileExtension = ".icseq";
        public static readonly string AnnotationFileExtension = ".icanno";
        public static readonly char Delimiter = '_';

        private static readonly ASCIIEncoding Encoding = new ASCIIEncoding();

        public FastaDatabase(string databaseFilePath)
        {
            if (!File.Exists(databaseFilePath))
                throw new FileNotFoundException("File not found: " + databaseFilePath);
            if (string.Compare(Path.GetExtension(databaseFilePath), ".fasta", StringComparison.OrdinalIgnoreCase) != 0)
                throw new FormatException("Not a fasta file: " + databaseFilePath);

            _databaseFilePath = databaseFilePath;
            _lastWriteTimeHash = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

            string databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);
            _seqFilePath = databaseFilePathNoExt + SeqFileExtension;
            _annoFilePath = databaseFilePathNoExt + AnnotationFileExtension;

            if (!File.Exists(_seqFilePath) || !File.Exists(_annoFilePath) || !CheckHashCode(_seqFilePath) || !CheckHashCode(_annoFilePath))
                GenerateMetaFiles();
        }

        public void Read()
        {
            if (!ReadSeqFile(_seqFilePath))
                throw new FormatException("Error while reading " + _seqFilePath);
            if (!ReadAnnnoFile(_annoFilePath))
                throw new FormatException("Error while reading " + _annoFilePath);
        }
            
        private readonly string _databaseFilePath;
        private readonly string _seqFilePath;
        private readonly string _annoFilePath;
        private readonly int _lastWriteTimeHash;
        private IDictionary<long, string> names;
        private IDictionary<long, string> descriptions;

        private void GenerateMetaFiles()
        {
            if(File.Exists(_seqFilePath))
                File.Delete(_seqFilePath);

            if(File.Exists(_annoFilePath))
                File.Delete(_annoFilePath);

            Console.WriteLine("Generating " + _seqFilePath + " and " + _annoFilePath + ".");

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
                    string name = reader.ProteinName;
                    string description = reader.ProteinDescription;
                    string sequence = Delimiter + reader.ProteinSequence;
                    seqWriter.Write(Encoding.GetBytes(sequence));
                    annoWriter.Write(offset + ":" + name + ":" + description);
                    offset += sequence.Length;
                }

                int hashCode = File.GetLastWriteTime(_databaseFilePath).GetHashCode();

                seqWriter.Write(hashCode);
                annoWriter.Write(hashCode);

                reader.CloseFile();
            }
        }

        private bool CheckHashCode(string filePath)
        {
            using (FileStream fs = File.OpenRead(filePath))
            {
                fs.Seek(-sizeof(int), SeekOrigin.End);

                using (var reader = new BinaryReader(fs))
                {
                    int lastWriteTimeHash = reader.ReadInt32();
                    if (lastWriteTimeHash == _lastWriteTimeHash)
                        return true;
                }
            }
            return false;
        }

        private bool ReadSeqFile(string seqFilePath)
        {
            using (var reader = new BinaryReader(File.Open(seqFilePath, FileMode.Open)))
            {
            }
            return true;
        }

        private bool ReadAnnnoFile(string annoFilePath)
        {
            names = new Dictionary<long, string>();
            descriptions = new Dictionary<long, string>();

            using (var reader = new StreamReader(_annoFilePath))
            {
                string s;
                while((s=reader.ReadLine()) != null)
                {
                    var token = s.Split(':');
                    if (token.Length != 3)
                        break;
                    var offset = long.Parse(token[0]);
                    names.Add(offset,token[1]);
                    descriptions.Add(offset, token[2]);
                }
            }

            return true;
        }
    }
}
