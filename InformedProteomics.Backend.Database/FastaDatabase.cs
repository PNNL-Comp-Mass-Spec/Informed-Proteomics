using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;
using ProteinFileReader;
using PRISM;

namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Facilitates working with a Fasta sequence database, using suffix arrays
    /// </summary>
    public class FastaDatabase
    {
        /// <summary>
        /// File format identifier to avoid incompatible backing files
        /// </summary>
        public const int FileFormatId = 175;

        /// <summary>
        /// Extension used for the Seq file
        /// </summary>
        public const string SeqFileExtension = ".icseq";

        /// <summary>
        /// Extension used for the Annotation file
        /// </summary>
        public const string AnnotationFileExtension = ".icanno";

        /// <summary>
        /// Extension used for the Decoy database file
        /// </summary>
        public const string DecoyDatabaseFileExtension = ".icdecoy.fasta";

        /// <summary>
        /// Extension used for the shuffled decoy database file
        /// </summary>
        public const string ShuffleDecoyFileExtension = ".icsfldecoy.fasta";

        /// <summary>
        /// For shuffled decoys, number of mutations
        /// </summary>
        public const int NumMutations = 3;

        /// <summary>
        /// Encoding used in the backing files
        /// </summary>
        public static readonly ASCIIEncoding Encoding = new ASCIIEncoding();

        /// <summary>
        /// True if this instance is tied to the decoy database
        /// </summary>
        public bool IsDecoy { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="databaseFilePath"></param>
        /// <param name="isDecoy"></param>
        public FastaDatabase(string databaseFilePath, bool isDecoy = false)
        {
            if (string.IsNullOrWhiteSpace(databaseFilePath))
            {
                throw new FileNotFoundException("Null or empty string passed to FastaDatabase");
            }

            if (!File.Exists(databaseFilePath))
            {
                throw new FileNotFoundException("File not found: " + databaseFilePath);
            }

            if (!FastaDatabaseConstants.ValidFASTAExtension(databaseFilePath))
            {
                throw new FormatException("Not a fasta file: " + databaseFilePath);
            }

            _databaseFilePath = databaseFilePath;
            _lastWriteTimeHash = File.GetLastWriteTimeUtc(_databaseFilePath).GetHashCode();

            _seqFilePath = Path.ChangeExtension(databaseFilePath, SeqFileExtension);
            _annoFilePath = Path.ChangeExtension(databaseFilePath, AnnotationFileExtension);

            if (!File.Exists(_seqFilePath)
                || !File.Exists(_annoFilePath)
                || !CheckHashCodeBinaryFile(_seqFilePath, _lastWriteTimeHash)
                || !CheckHashCodeTextFile(_annoFilePath, _lastWriteTimeHash))
            {
                Console.WriteLine("Generating " + _seqFilePath + " and ");
                Console.Write("Generating " + _annoFilePath + " ... ");

                DeleteOldIndexFiles(_databaseFilePath, isDecoy);

                GenerateMetaFiles();
                Console.WriteLine("Done");
            }

            IsDecoy = isDecoy;
        }

        /// <summary>
        /// Get the Decoy version of this database (create it if missing)
        /// </summary>
        /// <param name="enzyme"></param>
        /// <param name="shuffle"></param>
        /// <returns></returns>
        public FastaDatabase Decoy(Enzyme enzyme, bool shuffle = false)
        {
            if (IsDecoy)
            {
                throw new InvalidOperationException("Already a decoy database");
            }

            var decoyDatabasePath = GetDecoyDatabasePath(enzyme, shuffle);
            if (!File.Exists(decoyDatabasePath))
            {
                CreateDecoyDatabase(enzyme, shuffle);
            }

            return new FastaDatabase(decoyDatabasePath, true);
        }

        /// <summary>
        /// Create the decoy version of this database
        /// </summary>
        /// <param name="enzyme"></param>
        /// <param name="shuffle"></param>
        public void CreateDecoyDatabase(Enzyme enzyme, bool shuffle)
        {
            var reader = new FastaFileReader();

            var fastaFile = new FileInfo(_databaseFilePath);

            if (!fastaFile.Exists)
            {
                throw new FileNotFoundException("FASTA file not found: " + fastaFile.FullName);
            }

            if (!reader.OpenFile(_databaseFilePath))
            {
                return;
            }

            var decoyDatabaseFileName = GetDecoyDatabasePath(enzyme, shuffle);

            Console.WriteLine("Creating " + decoyDatabaseFileName);

            // Base the random number generator seed on the FASTA file size.
            // This assures that the same input FASTA file always gets the same shuffled sequences
            var seed = fastaFile.Length;

            while (seed > int.MaxValue)
            {
                seed = (long)Math.Floor(seed / 2.0);
            }

            SimpleStringProcessing.DefineRandomNumberGeneratorSeed((int)seed);

            using (var decoyWriter = new StreamWriter(decoyDatabaseFileName))
            {
                while (reader.ReadNextProteinEntry())
                {
                    var name = reader.ProteinName;
                    var description = reader.ProteinDescription;
                    var sequence = reader.ProteinSequence;

                    decoyWriter.WriteLine(">{0}_{1} {2}", FastaDatabaseConstants.DecoyProteinPrefix, name, description);

                    if (!shuffle)
                    {
                        // Reversed protein sequence
                        var decoySequence = new StringBuilder();
                        for (var i = sequence.Length - 1; i >= 0; i--)
                        {
                            var residue = sequence[i];
                            if (enzyme?.Residues.Length > 0 && enzyme.IsCleavable(residue) && decoySequence.Length > 0)
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
                        // Shuffled protein sequences
                        decoyWriter.WriteLine(SimpleStringProcessing.Mutate(SimpleStringProcessing.Shuffle(sequence), NumMutations));
                    }
                }
                reader.CloseFile();
            }
        }

        /// <summary>
        /// Returns the characters in the sequence file
        /// </summary>
        /// <returns></returns>
        public IEnumerable<byte> Characters()
        {
            if (_sequence != null)
            {
                for (var i = 0; i < _sequence.Length - 1; i++)
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
        }

        //public IEnumerable<byte> Characters(int index, int numCharacters)
        //{
        //    if (_sequence == null) throw new SystemException("Fasta sequence must be loaded!");
        //    for (var i = index; i < index+numCharacters && i < _sequence.Length - 1; i++)
        //    {
        //        yield return _sequence[i];
        //    }
        //}

        /// <summary>
        /// Read in the backing files
        /// </summary>
        public void Read()
        {
            if (!ReadSeqFile())
            {
                throw new FormatException("Error while reading " + _seqFilePath);
            }

            if (!ReadAnnoFile())
            {
                throw new FormatException("Error while reading " + _annoFilePath);
            }
        }

        /// <summary>
        /// Path to the Fasta file
        /// </summary>
        /// <returns></returns>
        public string GetFastaFilePath()
        {
            return _databaseFilePath;
        }

        /// <summary>
        /// Number of proteins in the database
        /// </summary>
        /// <returns></returns>
        public int GetNumEntries()
        {
            Read();
            return _names.Count;
        }

        /// <summary>
        /// Get the names of all proteins in the database
        /// </summary>
        /// <returns></returns>
        public IEnumerable<string> GetProteinNames()
        {
            Read();
            return _names.Values;
        }

        /// <summary>
        /// Get the entire concatenated sequence
        /// </summary>
        /// <returns></returns>
        public byte[] GetSequence()
        {
            if (_sequence == null)
            {
                Read();
            }

            return _sequence;
        }

        /// <summary>
        /// Print the entire sequence to console
        /// </summary>
        public void PrintSequence()
        {
            Console.WriteLine(_sequence == null ? "Annotation is null!" : Encoding.GetString(_sequence));
        }

        /// <summary>
        /// Generate the path for the decoy database according to the supplied parameters
        /// </summary>
        /// <param name="enzyme"></param>
        /// <param name="shuffle"></param>
        /// <returns></returns>
        public string GetDecoyDatabasePath(Enzyme enzyme, bool shuffle = false)
        {
            string newExtension;

            if (!shuffle)
            {
                // Reverse the sequences
                if (enzyme?.Residues.Length > 0)
                {
                    newExtension = ".icdecoy." + new string(enzyme.Residues) + ".fasta";
                }
                else
                {
                    newExtension = DecoyDatabaseFileExtension;
                }
            }
            else
            {
                // Shuffle the sequences
                newExtension = ShuffleDecoyFileExtension;
            }

            return Path.ChangeExtension(_databaseFilePath, newExtension);
        }

        /// <summary>
        /// Get the name of the protein that starts at <paramref name="offset"/> in the concatenated sequence
        /// </summary>
        /// <param name="offset"></param>
        /// <returns></returns>
        public string GetProteinName(long offset)
        {
            var offsetKey = GetOffsetKey(offset);
            if (_names.TryGetValue(offsetKey, out var proteinName))
            {
                return proteinName;
            }

            return "UnknownProtein_Offset" + offset;
        }

        /// <summary>
        /// Get the description of the protein that starts at <paramref name="offset"/> in the concatenated sequence
        /// </summary>
        /// <param name="offset"></param>
        /// <returns></returns>
        public string GetProteinDescription(long offset)
        {
            var offsetKey = GetOffsetKey(offset);

            if (_descriptions.TryGetValue(offsetKey, out var proteinDescription))
            {
                return proteinDescription;
            }

            return "Unknown description, Offset " + offset;
        }

        /// <summary>
        /// Get the offset in the concatenated sequence of the protein with name <paramref name="name"/>
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public long? GetOffset(string name)
        {
            if (!_nameToOffset.TryGetValue(name, out var offset))
            {
                return null;
            }

            return offset;
        }

        /// <summary>
        /// Get the description of the protein with name <paramref name="name"/>
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public string GetProteinDescription(string name)
        {
            if (!_nameToOffset.TryGetValue(name, out var offset))
            {
                return null;
            }

            var offsetKey = GetOffsetKey(offset);
            return _descriptions[offsetKey];
        }

        /// <summary>
        /// Get the length of the protein with name <paramref name="name"/>
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public int GetProteinLength(string name)
        {
            if (_nameToLength.TryGetValue(name, out var length))
            {
                return length;
            }

            return -1;
        }

        /// <summary>
        /// Get the sequence of the protein with name <paramref name="name"/>
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public string GetProteinSequence(string name)
        {
            if (!_nameToOffset.TryGetValue(name, out var offset))
            {
                return null;
            }

            var length = _nameToLength[name];
            return GetProteinSequence(offset, length);
        }

        /// <summary>
        /// Return the protein sequence starting at the given offset, spanning the given length
        /// </summary>
        /// <param name="offset"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private string GetProteinSequence(long offset, int length)
        {
            return Encoding.GetString(_sequence, (int)(offset + 1), length);
        }

        /// <summary>
        /// Get the position in the protein sequence of the offset in the concatenated sequence (one-based index)
        /// </summary>
        /// <param name="offset"></param>
        /// <returns></returns>
        public int GetOneBasedPositionInProtein(long offset)
        {
            return (int)(offset - GetOffsetKey(offset));
        }

        /// <summary>
        /// Get the position in the protein sequence of the offset in the concatenated sequence (zero-based index)
        /// </summary>
        /// <param name="offset"></param>
        /// <returns></returns>
        public int GetZeroBasedPositionInProtein(long offset)
        {
            return GetOneBasedPositionInProtein(offset) - 1;
        }

        /// <summary>
        /// Returns the hash based on last write time, used for consistency verification
        /// </summary>
        /// <returns></returns>
        internal int GetLastWriteTimeHash()
        {
            return _lastWriteTimeHash;
        }

        /// <summary>
        /// For file <paramref name="filePath"/>, check the last write time hash against <paramref name="code"/>
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="code"></param>
        /// <returns></returns>
        internal static bool CheckHashCodeBinaryFile(string filePath, int code)
        {
            var dataFile = new FileInfo(filePath);
            if (dataFile.Length < 2 * sizeof(int))
            {
                return false;
            }

            using (var fs = new FileStream(dataFile.FullName, FileMode.Open, FileAccess.Read, FileShare.Read))
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-2 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId != FileFormatId)
                {
                    return false;
                }

                var lastWriteTimeHash = reader.ReadInt32();
                if (lastWriteTimeHash == code)
                {
                    return true;
                }
            }

            return false;
        }

        internal static bool CheckHashCodeTextFile(string filePath, int code)
        {
            var dataFile = new FileInfo(filePath);
            if (dataFile.Length < 1)
            {
                return false;
            }

            var lastLine = string.Empty;

            using (var fs = new FileStream(dataFile.FullName, FileMode.Open, FileAccess.Read, FileShare.Read))
            using (var reader = new StreamReader(fs))
            {
                while (!reader.EndOfStream)
                {
                    lastLine = reader.ReadLine();
                }
            }

            if (string.IsNullOrWhiteSpace(lastLine))
            {
                return false;
            }

            var token = lastLine.Split(FastaDatabaseConstants.AnnotationDelimiter);
            if (token.Length != 2)
            {
                return false;
            }

            var fileFormatId = Convert.ToInt32(token[0]);
            if (fileFormatId != FileFormatId)
            {
                return false;
            }

            var lastWriteTimeHash = Convert.ToInt32(token[1]);
            if (lastWriteTimeHash == code)
            {
                return true;
            }

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

        /// <summary>
        /// Tracks duplicate names along with a count for each name
        /// </summary>
        /// <remarks>Used to auto-rename proteins</remarks>
        private IDictionary<string, int> _duplicateNameCounts;

        private byte[] _sequence;

        private void DeleteOldIndexFiles(string dbFilePath, bool isDecoy)
        {
            var dbFile = new FileInfo(dbFilePath);
            var dataset = Path.GetFileNameWithoutExtension(dbFile.Name);

            var diFolder = dbFile.Directory;
            if (diFolder == null)
            {
                return;
            }

            var fileNamesToDelete = new List<string>
            {
                dataset + ".icdecoy*",
                // ReSharper disable StringLiteralTypo
                dataset + ".icsfldecoy*"
                // ReSharper restore StringLiteralTypo
            };

            if (!isDecoy)
            {
                fileNamesToDelete.Add(dataset + SeqFileExtension);            // .icseq
                fileNamesToDelete.Add(dataset + AnnotationFileExtension);      // .icanno
                fileNamesToDelete.Add(dataset + ".icplcp");
            }

            foreach (var fileName in fileNamesToDelete)
            {
                foreach (var fileToDelete in diFolder.GetFiles(fileName))
                {
                    try
                    {
                        fileToDelete.Delete();
                    }
                    catch (Exception ex)
                    {
                        throw new Exception("Error deleting index file " + fileToDelete.Name, ex);
                    }
                }
            }
        }

        private void GenerateMetaFiles()
        {
            if (File.Exists(_seqFilePath))
            {
                File.Delete(_seqFilePath);
            }

            if (File.Exists(_annoFilePath))
            {
                File.Delete(_annoFilePath);
            }

            // Keys are protein name
            // Values track the number of times the name has been encountered, the number of residues, and a SHA1 hash of the residues
            var proteinNamesAndStats = new Dictionary<string, ProteinHashInfo>(StringComparer.InvariantCultureIgnoreCase);

            using (var seqWriter = new BinaryWriter(File.Open(_seqFilePath, FileMode.CreateNew)))
            using (var annoWriter = new StreamWriter(_annoFilePath))
            {
                // Read
                var reader = new FastaFileReader();
                if (!reader.OpenFile(_databaseFilePath))
                {
                    return;
                }

                long offset = 0;
                while (reader.ReadNextProteinEntry())
                {
                    var name = reader.ProteinName;
                    var description = reader.ProteinDescription;
                    var sequence = (char) FastaDatabaseConstants.Delimiter + reader.ProteinSequence;
                    var length = sequence.Length;

                    var proteinInfoCurrent = new ProteinHashInfo(sequence);

                    if (proteinNamesAndStats.TryGetValue(name, out var proteinInfoCached))
                    {
                        // Duplicate name; either skip this protein or rename it
                        if (proteinInfoCached.SequenceLength == proteinInfoCurrent.SequenceLength &&
                            string.Equals(proteinInfoCached.SequenceHash, proteinInfoCurrent.SequenceHash))
                        {
                            // Duplicate protein; skip it
                            continue;
                        }

                        name += "_Duplicate" + proteinInfoCached.ObservationCount.ToString("00");
                        proteinInfoCached.ObservationCount++;
                    }

                    proteinNamesAndStats.Add(name, proteinInfoCurrent);

                    seqWriter.Write(Encoding.GetBytes(sequence));

                    annoWriter.WriteLine("{0}{1}{2}{3}{4}{5}{6}",
                        offset, FastaDatabaseConstants.AnnotationDelimiter,
                        length, FastaDatabaseConstants.AnnotationDelimiter,
                        name.Replace(FastaDatabaseConstants.AnnotationDelimiter, '-'), FastaDatabaseConstants.AnnotationDelimiter,
                        description);

                    offset += sequence.Length;
                }

                seqWriter.Write((char) FastaDatabaseConstants.Delimiter);

                // write file format Id
                seqWriter.Write(FileFormatId);

                // writer last write time hash
                var hashCode = File.GetLastWriteTimeUtc(_databaseFilePath).GetHashCode();

                seqWriter.Write(hashCode);
                annoWriter.Write("{0}{1}{2}", FileFormatId, FastaDatabaseConstants.AnnotationDelimiter, hashCode);

                reader.CloseFile();
            }
        }

        private bool ReadSeqFile()
        {
            using (var fileStream = new FileStream(_seqFilePath, FileMode.Open, FileAccess.Read))
            {
                _sequence = new byte[fileStream.Length - 2 * sizeof(int) + 1];
                fileStream.Read(_sequence, 0, _sequence.Length);
                _sequence[_sequence.Length - 1] = FastaDatabaseConstants.LastCharacter;
            }

            return true;
        }

        private bool ReadAnnoFile()
        {
            _offsetList = new List<long>();
            _nameToLength = new Dictionary<string, int>();
            _names = new Dictionary<long, string>();
            _descriptions = new Dictionary<long, string>();
            _nameToOffset = new Dictionary<string, long>();

            _duplicateNameCounts= new Dictionary<string, int>();

            using (var reader = new StreamReader(_annoFilePath))
            {
                string s;
                while ((s = reader.ReadLine()) != null)
                {
                    var token = s.Split(FastaDatabaseConstants.AnnotationDelimiter);
                    if (token.Length < 4)
                    {
                        break;
                    }

                    var offset = long.Parse(token[0]);
                    _offsetList.Add(offset);
                    var length = int.Parse(token[1]);
                    var name = token[2];

                    if (_nameToLength.TryGetValue(name, out _))
                    {
                        ConsoleMsgUtils.ShowWarning(string.Format(
                            "Duplicate protein name, renaming {0} at offset {1} in {2} to avoid collisions",
                            name, offset, Path.GetFileName(_annoFilePath)));

                        if (_duplicateNameCounts.TryGetValue(name, out var duplicateSuffix))
                        {
                            duplicateSuffix++;
                            _duplicateNameCounts[name] = duplicateSuffix;
                        }
                        else
                        {
                            duplicateSuffix = 1;
                            _duplicateNameCounts.Add(name, duplicateSuffix);
                        }
                        name += "_Duplicate" + duplicateSuffix.ToString("00");
                    }

                    _nameToLength.Add(name, length);
                    _names.Add(offset, name);
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
