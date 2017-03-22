using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace MSPathFinderT
{
    using InformedProteomics.Backend.Data.Spectrometry;

    public class TopDownInputParameters
    {
        public const string ParameterFileExtension = ".param";

        public IEnumerable<string> SpecFilePaths { get; set; }
        public string DatabaseFilePath { get; set; }
        public string OutputDir { get; set; }
        public AminoAcidSet AminoAcidSet { get; set; }

        public int SearchModeInt
        {
            get
            {
                if (SearchMode == InternalCleavageType.MultipleInternalCleavages)
                    return 0;
                if (SearchMode == InternalCleavageType.SingleInternalCleavage)
                    return 1;
                return 2;
            }
            set
            {
                if (value == 0)
                {
                    SearchMode = InternalCleavageType.MultipleInternalCleavages;
                }
                else if (value == 1)
                {
                    SearchMode = InternalCleavageType.SingleInternalCleavage;
                }
                else
                {
                    SearchMode = InternalCleavageType.NoInternalCleavage;
                }
            }
        }

        public InternalCleavageType SearchMode { get; set; }

        public bool TagBasedSearch { get; set; }

        public bool? TdaBool
        {
            get
            {
                if (Tda == DatabaseSearchMode.Both)
                    return true;
                if (Tda == DatabaseSearchMode.Decoy)
                    return null;
                //(Tda2 == DatabaseSearchMode.Target)
                return false;
            }
            set
            {
                if (value == null)
                {
                    Tda = DatabaseSearchMode.Decoy;
                }
                else if (value.Value)
                {
                    Tda = DatabaseSearchMode.Both;
                }
                else
                {
                    Tda = DatabaseSearchMode.Target;
                }
            }
        }

        public DatabaseSearchMode Tda { get; set; }
        public double PrecursorIonTolerancePpm { get; set; }
        public double ProductIonTolerancePpm { get; set; }
        public int MinSequenceLength { get; set; }
        public int MaxSequenceLength { get; set; }
        public int MinPrecursorIonCharge { get; set; }
        public int MaxPrecursorIonCharge { get; set; }
        public int MinProductIonCharge { get; set; }
        public int MaxProductIonCharge { get; set; }
        public double MinSequenceMass { get; set; }
        public double MaxSequenceMass { get; set; }
        public string FeatureFilePath { get; set; }
        public int MaxNumThreads { get; set; }
        public ActivationMethod ActivationMethod { get; set; }

        /// <summary>
        /// Specific MS2 scans to process
        /// </summary>
        public SortedSet<int> ScanNumbers { get; set; }

        private IEnumerable<SearchModification> _searchModifications;
        private int _maxNumDynModsPerSequence;

        public void Display()
        {
            Console.WriteLine("MaxThreads: " + MaxNumThreads);

            foreach (var specFilePath in SpecFilePaths)
            {
                Console.WriteLine("SpectrumFilePath: " + specFilePath);
            }

            Console.WriteLine("DatabaseFilePath: " + DatabaseFilePath);
            Console.WriteLine("FeatureFilePath:  {0}", FeatureFilePath ?? "N/A");
            Console.WriteLine("OutputDir:        " + OutputDir);
            Console.WriteLine("SearchMode: " + SearchModeInt);
            Console.WriteLine("Tag-based search: " + TagBasedSearch);
            Console.WriteLine("Tda: " + (TdaBool == null ? "Decoy" : (bool)TdaBool ? "Target+Decoy" : "Target"));
            Console.WriteLine("PrecursorIonTolerancePpm: " + PrecursorIonTolerancePpm);
            Console.WriteLine("ProductIonTolerancePpm: " + ProductIonTolerancePpm);
            Console.WriteLine("MinSequenceLength: " + MinSequenceLength);
            Console.WriteLine("MaxSequenceLength: " + MaxSequenceLength);
            Console.WriteLine("MinPrecursorIonCharge: " + MinPrecursorIonCharge);
            Console.WriteLine("MaxPrecursorIonCharge: " + MaxPrecursorIonCharge);
            Console.WriteLine("MinProductIonCharge: " + MinProductIonCharge);
            Console.WriteLine("MaxProductIonCharge: " + MaxProductIonCharge);
            Console.WriteLine("MinSequenceMass: " + MinSequenceMass);
            Console.WriteLine("MaxSequenceMass: " + MaxSequenceMass);
            Console.WriteLine("MaxDynamicModificationsPerSequence: " + _maxNumDynModsPerSequence);
            Console.WriteLine("Modifications: ");

            foreach (var searchMod in _searchModifications)
            {
                Console.WriteLine(searchMod);
            }

            if (!string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                Console.WriteLine("Getting MS1 features from " + FeatureFilePath);
            }

            if (ScanNumbers != null && ScanNumbers.Any())
            {
                Console.WriteLine("Processing specific MS2 scans:");
                Console.WriteLine(string.Join(", ", ScanNumbers));
            }
        }

        public void Write()
        {
            foreach (var specFilePath in SpecFilePaths)
            {
                var outputFilePath = Path.Combine(OutputDir, Path.GetFileNameWithoutExtension(specFilePath) + ParameterFileExtension);

                using (var writer = new StreamWriter(outputFilePath))
                {
                    writer.WriteLine("SpecFile\t" + Path.GetFileName(specFilePath));
                    writer.WriteLine("DatabaseFile\t" + Path.GetFileName(DatabaseFilePath));
                    writer.WriteLine("FeatureFile\t{0}", FeatureFilePath != null ? Path.GetFileName(FeatureFilePath) : Path.GetFileName(MassSpecDataReaderFactory.ChangeExtension(specFilePath, ".ms1ft")));
                    writer.WriteLine("SearchMode\t" + SearchModeInt);
                    writer.WriteLine("Tag-based search\t" + TagBasedSearch);
                    writer.WriteLine("Tda\t" + (TdaBool == null ? "Decoy" : (bool)TdaBool ? "Target+Decoy" : "Target"));
                    writer.WriteLine("PrecursorIonTolerancePpm\t" + PrecursorIonTolerancePpm);
                    writer.WriteLine("ProductIonTolerancePpm\t" + ProductIonTolerancePpm);
                    writer.WriteLine("MinSequenceLength\t" + MinSequenceLength);
                    writer.WriteLine("MaxSequenceLength\t" + MaxSequenceLength);
                    writer.WriteLine("MinPrecursorIonCharge\t" + MinPrecursorIonCharge);
                    writer.WriteLine("MaxPrecursorIonCharge\t" + MaxPrecursorIonCharge);
                    writer.WriteLine("MinProductIonCharge\t" + MinProductIonCharge);
                    writer.WriteLine("MaxProductIonCharge\t" + MaxProductIonCharge);
                    writer.WriteLine("MinSequenceMass\t" + MinSequenceMass);
                    writer.WriteLine("MaxSequenceMass\t" + MaxSequenceMass);
                    writer.WriteLine("MaxDynamicModificationsPerSequence\t" + _maxNumDynModsPerSequence);
                    foreach (var searchMod in _searchModifications)
                    {
                        writer.WriteLine("Modification\t" + searchMod);
                    }
                }
            }
        }

        public string Parse(Dictionary<string, string> parameters)
        {
            var message = CheckIsValid(parameters);
            if (message != null)
                return message;

            try
            {
                var specFilePath = parameters["-s"];
                SpecFilePaths = Directory.Exists(specFilePath) ? Directory.GetFiles(specFilePath, "*.raw") : new[] { specFilePath };
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception parsing the file path for parameter -s: " + ex.Message);
                throw;
            }

            DatabaseFilePath = parameters["-d"];

            try
            {
                var outputDir = parameters["-o"] ?? Environment.CurrentDirectory;

                if (!Directory.Exists(outputDir))
                {
                    if (File.Exists(outputDir) && !File.GetAttributes(outputDir).HasFlag(FileAttributes.Directory))
                    {
                        return "OutputDir is not a directory: " + outputDir;
                    }
                    Directory.CreateDirectory(outputDir);
                }
                OutputDir = outputDir;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception validating the path for parameter -o: " + ex.Message);
                throw;
            }

            try
            {
                var modFilePath = parameters["-mod"];
                var errorMessage = LoadModsFile(modFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                    return errorMessage;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception parsing the file for parameter -mod: " + ex.Message);
                throw;
            }

            try
            {
                var scansFilePath = parameters["-scansFile"];
                var errorMessage = LoadScansFile(scansFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                    return errorMessage;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception parsing the file for parameter -scansFile: " + ex.Message);
                throw;
            }

            var currentParameter = string.Empty;

            try
            {
                currentParameter = "-feature";
                FeatureFilePath = parameters["-feature"];

                currentParameter = "-m";
                var searchMode = Convert.ToInt32(parameters["-m"]);
                if (searchMode < 0 || searchMode > 2)
                {
                    return "Invalid value (" + searchMode + ") for parameter -m";
                }
                SearchModeInt = searchMode;

                currentParameter = "-t";
                PrecursorIonTolerancePpm = Convert.ToDouble(parameters["-t"]);

                currentParameter = "-f";
                ProductIonTolerancePpm = Convert.ToDouble(parameters["-f"]);

                currentParameter = "-tda";
                var tdaVal = Convert.ToInt32(parameters["-tda"]);
                if (tdaVal != 0 && tdaVal != 1 && tdaVal != -1)
                {
                    return "Invalid value (" + tdaVal + ") for parameter -tda";
                }

                if (tdaVal == 1)
                {
                    Tda = DatabaseSearchMode.Both;
                }
                else if (tdaVal == -1)
                {
                    Tda = DatabaseSearchMode.Decoy;
                }
                else
                {
                    Tda = DatabaseSearchMode.Target;
                }

                currentParameter = "-minLength";
                MinSequenceLength = Convert.ToInt32(parameters["-minLength"]);

                currentParameter = "-maxLength";
                MaxSequenceLength = Convert.ToInt32(parameters["-maxLength"]);
                if (MinSequenceLength > MaxSequenceLength)
                {
                    return "MinSequenceLength (" + MinSequenceLength + ") is larger than MaxSequenceLength (" + MaxSequenceLength + ")!";
                }

                currentParameter = "-feature";
                MinPrecursorIonCharge = Convert.ToInt32(parameters["-minCharge"]);

                currentParameter = "-feature";
                MaxPrecursorIonCharge = Convert.ToInt32(parameters["-maxCharge"]);
                if (MinSequenceLength > MaxSequenceLength)
                {
                    return "MinPrecursorCharge (" + MinPrecursorIonCharge + ") is larger than MaxPrecursorCharge (" + MaxPrecursorIonCharge + ")!";
                }

                currentParameter = "-minFragCharge";
                MinProductIonCharge = Convert.ToInt32(parameters["-minFragCharge"]);

                currentParameter = "-maxFragCharge";
                MaxProductIonCharge = Convert.ToInt32(parameters["-maxFragCharge"]);
                if (MinSequenceLength > MaxSequenceLength)
                {
                    return "MinFragmentCharge (" + MinProductIonCharge + ") is larger than MaxFragmentCharge (" + MaxProductIonCharge + ")!";
                }

                currentParameter = "-minMass";
                MinSequenceMass = Convert.ToDouble(parameters["-minMass"]);

                currentParameter = "-maxMass";
                MaxSequenceMass = Convert.ToDouble(parameters["-maxMass"]);
                if (MinSequenceMass > MaxSequenceMass)
                {
                    return "MinSequenceMassInDa (" + MinSequenceMass + ") is larger than MaxSequenceMassInDa (" + MaxSequenceMass + ")!";
                }

                currentParameter = "-threads";
                MaxNumThreads = Convert.ToInt32(parameters["-threads"]);
                MaxNumThreads = GetOptimalMaxThreads(MaxNumThreads);

                currentParameter = "-tagSearch";
                TagBasedSearch = Convert.ToInt32(parameters["-tagSearch"]) == 1;

                currentParameter = "-act";
                ActivationMethod = (ActivationMethod)Convert.ToInt32(parameters["-act"]);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception parsing parameter '" + currentParameter + "': " + ex.Message);
                throw;
            }

            return null;
        }

        static string CheckIsValid(Dictionary<string, string> parameters)
        {
            var optionalKeys = new SortedSet<string>(StringComparer.InvariantCultureIgnoreCase)
            {
                "-mod",
                "-o",
                "-feature",
                "-scansFile"
            };

            foreach (var keyValuePair in parameters)
            {
                var key = keyValuePair.Key;
                var value = keyValuePair.Value;

                try
                {
                    if (string.IsNullOrWhiteSpace(keyValuePair.Value) && !optionalKeys.Contains(keyValuePair.Key))
                    {
                        return "Missing required parameter " + key;
                    }

                    if (key.Equals("-s"))
                    {
                        if (value == null)
                        {
                            return "Missing parameter " + key;
                        }
                        if (!File.Exists(value) && !Directory.Exists(value))
                        {
                            return "File not found: " + value;
                        }
                        if (Directory.Exists(value)) continue;
                        if (MassSpecDataReaderFactory.GetMassSpecDataReader(value) == null)
                        {
                            var extension = Path.GetExtension(value);
                            return "Invalid extension for the parameter " + key + " (" + extension + ")";
                        }
                    }
                    else if (key.Equals("-d"))
                    {
                        if (value == null)
                        {
                            return "Missing required parameter " + key;
                        }
                        if (!File.Exists(value))
                        {
                            return "File not found: " + value;
                        }
                        var extension = Path.GetExtension(value).ToLower();
                        if (!extension.Equals(".fa") &&
                            !extension.Equals(".fasta"))
                        {
                            return "Invalid extension for the parameter " + key + " (" + extension + ")";
                        }
                    }
                    else if (key.Equals("-o"))
                    {
                    }
                    else if (key.Equals("-threads"))
                    {
                        if (value == null || Int32.Parse(value) < 0)
                        {
                            return "Invalid number of maximum threads " + value;
                        }
                    }
                    else if (key.Equals("-mod"))
                    {
                        if (!string.IsNullOrWhiteSpace(value) && !File.Exists(value))
                        {
                            return "-mod file not found: " + value;
                        }
                    }
                    else if (key.Equals("-scansFile"))
                    {
                        if (!string.IsNullOrWhiteSpace(value) && !File.Exists(value))
                        {
                            return "-scansFile file not found: " + value;
                        }
                    }
                    else if (key.Equals("-feature"))
                    {
                        if (value != null && !File.Exists(value))
                        {
                            return "File not found: " + value;
                        }
                        if (value != null &&
                            !Path.GetExtension(value).ToLower().Equals(".csv") &&
                            !Path.GetExtension(value).ToLower().Equals(".ms1ft") &&
                            !Path.GetExtension(value).ToLower().Equals(".msalign"))
                        {
                            return "Invalid extension for the parameter " + key + " (" + Path.GetExtension(value) + ")";
                        }
                    }
                    else if (key.Equals("-act"))
                    {
                        int val;
                        if (!Int32.TryParse(value, out val) || val < 0 || val > 6)
                        {
                            return "Invalid value for the parameter " + key + " {" + value + ")";
                        }
                    }
                    else
                    {
                        double num;
                        if (!double.TryParse(value, out num))
                        {
                            return "Invalid value (" + value + ") for the parameter " + key;
                        }
                    }
                }
                catch (Exception)
                {
                    if (string.IsNullOrEmpty(key))
                        key = "?UnknownKey?";

                    if (string.IsNullOrEmpty(value))
                        value = string.Empty;

                    Console.WriteLine("Error parsing parameter '" + key + "' with value '" + value + "'");
                    throw;
                }
            }

            return null;
        }

        /// <summary>
        /// Prevent the user from specifying too many threads
        /// </summary>
        /// <param name="userMaxThreads"></param>
        /// <returns></returns>
        private int GetOptimalMaxThreads(int userMaxThreads)
        {
            var threads = userMaxThreads;

            // If user-supplied thread count is zero or greater than NumLogicalCores, set it dynamically
            if (threads <= 0 || threads > ParallelizationUtils.NumLogicalCores)
            {
                // Non-hyperthreaded max: NumPhysicalCores, which will be the same as NumLogicalCores
                threads = ParallelizationUtils.NumPhysicalCores;

                // If hyperthreaded: Don't max out the threads, drop it down by one or two to maintain system responsiveness
                if (ParallelizationUtils.NumLogicalCores > ParallelizationUtils.NumPhysicalCores)
                {
                    if (ParallelizationUtils.NumLogicalCores - 2 > ParallelizationUtils.NumPhysicalCores)
                    {
                        threads = ParallelizationUtils.NumLogicalCores - 2;
                    }
                    else
                    {
                        threads = ParallelizationUtils.NumLogicalCores - 1;
                    }
                }
            }
            return threads;
        }

        private string LoadModsFile(string modFilePath)
        {
            if (string.IsNullOrWhiteSpace(modFilePath))
            {
                AminoAcidSet = new AminoAcidSet();
                _searchModifications = new SearchModification[0];
                return string.Empty;
            }

            if (!File.Exists(modFilePath))
            {
                return "-mod file not found: " + modFilePath;
            }

            var parser = new ModFileParser(modFilePath);
            _searchModifications = parser.SearchModifications;
            _maxNumDynModsPerSequence = parser.MaxNumDynModsPerSequence;

            if (_searchModifications == null)
                return "Error while parsing " + modFilePath;

            AminoAcidSet = new AminoAcidSet(_searchModifications, _maxNumDynModsPerSequence);
            return string.Empty;
        }

        private string LoadScansFile(string scansFilePath)
        {
            if (string.IsNullOrWhiteSpace(scansFilePath))
                return string.Empty;

            if (!File.Exists(scansFilePath))
            {
                return "-scansFile not found: " + scansFilePath;
            }

            ScanNumbers = new SortedSet<int>();

            var delimiters = new[] { ' ', '\t', ',' };

            using (var reader = new StreamReader(new FileStream(scansFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                while (!reader.EndOfStream)
                {
                    var dataLine = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(dataLine))
                        continue;

                    var dataValues = dataLine.Split(delimiters);
                    foreach (var value in dataValues)
                    {
                        int scanNumber;
                        if (!int.TryParse(value, out scanNumber))
                            continue;

                        if (!ScanNumbers.Contains(scanNumber))
                            ScanNumbers.Add(scanNumber);
                    }
                }
            }

            return string.Empty;
        }
    }
}
