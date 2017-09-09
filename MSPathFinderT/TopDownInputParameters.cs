using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using PRISM;

namespace MSPathFinderT
{
    public class TopDownInputParameters
    {
        public const string ParameterFileExtension = ".param";

        [Option("s", "specFile", Required = true, HelpText = "Spectrum File (*.raw)", HelpShowsDefault = false)]
        public string SpecFilePath { get; set; }

        public IEnumerable<string> SpecFilePaths { get; set; }

        [Option("d", "database", Required = true, HelpText = "Database File (*.fasta or *.fa)", HelpShowsDefault = false)]
        public string DatabaseFilePath { get; set; }

        [Option("o", "outputDir", HelpText = "Output Folder", HelpShowsDefault = false)]
        public string OutputDir { get; set; }

        public AminoAcidSet AminoAcidSet { get; set; }

        [Option("m", "searchMode", Min = 0, Max = 2, HelpText = "Search Mode (old format) (0: multiple internal cleavages, 1: single internal cleavage, 2: no internal cleavage)")]
        public int SearchModeInt
        {
            get
            {
                if (InternalCleavageMode == InternalCleavageType.MultipleInternalCleavages)
                    return 0;
                if (InternalCleavageMode == InternalCleavageType.SingleInternalCleavage)
                    return 1;
                return 2;
            }
            set
            {
                if (value == 0)
                {
                    InternalCleavageMode = InternalCleavageType.MultipleInternalCleavages;
                }
                else if (value == 1)
                {
                    InternalCleavageMode = InternalCleavageType.SingleInternalCleavage;
                }
                else
                {
                    InternalCleavageMode = InternalCleavageType.NoInternalCleavage;
                }
            }
        }

        [Option("ic", HelpText = "Search Mode")]
        public InternalCleavageType InternalCleavageMode { get; set; }

        [Option("tagSearch", HelpText = "Include Tag-based Search (can use '0' for false, '1' for true)")]
        public bool TagBasedSearch { get; set; }

        [Option("mod", HelpText = "Path to modification file. (Default: no modifications)", HelpShowsDefault = false)]
        public string ModsFilePath { get; set; }

        [Option("tda", Min = -1, Max = 1, HelpText = "Database search mode (0: don't search decoy database, 1: search shuffled decoy database)")]
        public int TdaInt
        {
            get
            {
                if (TargetDecoySearchMode == DatabaseSearchMode.Both)
                    return 1;
                if (TargetDecoySearchMode == DatabaseSearchMode.Decoy)
                    return -1;
                //(Tda2 == DatabaseSearchMode.Target)
                return 0;
            }
            set
            {
                if (value == -1)
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Decoy;
                }
                else if (value == 1)
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Both;
                }
                else
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Target;
                }
            }
        }

        //[Option("tda", HelpText = "Database search mode")]
        public DatabaseSearchMode TargetDecoySearchMode { get; set; }

        [Option("t", "precursorTol", /*Min = 1,*/ HelpText = "Precursor Tolerance (in PPM)")]
        public double PrecursorIonTolerancePpm { get; set; }

        [Option("f", "fragmentTol", /*Min = 1,*/ HelpText = "Fragment Ion Tolerance (in PPM)")]
        public double ProductIonTolerancePpm { get; set; }

        [Option("minLength", Min = 0, HelpText = "Minimum Sequence Length")]
        public int MinSequenceLength { get; set; }

        [Option("maxLength", Min = 0, HelpText = "Maximum Sequence Length")]
        public int MaxSequenceLength { get; set; }

        [Option("minCharge", Min = 1, HelpText = "Minimum precursor ion charge")]
        public int MinPrecursorIonCharge { get; set; }

        [Option("maxCharge", Min = 1, HelpText = "Maximum precursor ion charge")]
        public int MaxPrecursorIonCharge { get; set; }

        [Option("minFragCharge", Min = 1, HelpText = "Minimum fragment ion charge")]
        public int MinProductIonCharge { get; set; }

        [Option("maxFragCharge", Min = 1, HelpText = "Maximum fragment ion charge")]
        public int MaxProductIonCharge { get; set; }

        [Option("minMass", /*Min = 1,*/ HelpText = "Minimum sequence mass in Da")]
        public double MinSequenceMass { get; set; }

        [Option("maxMass", /*Min = 1,*/ HelpText = "Maximum sequence mass in Da")]
        public double MaxSequenceMass { get; set; }

        [Option("feature", HelpText = "*.ms1ft, *_isos.csv, or *.msalign, (Default: Run ProMex)", HelpShowsDefault = false)]
        public string FeatureFilePath { get; set; }

        [Option("threads", Min = 0, HelpText = "Maximum number of threads, or 0 to set automatically")]
        public int MaxNumThreads { get; set; }

        [Option("act", HelpText = "Activation Method")]
        public ActivationMethod ActivationMethod { get; set; }

        [Option("scansFile", HelpText = "Text file with MS2 scans to process", HelpShowsDefault = false)]
        public string ScansFilePath { get; set; }

        public TopDownInputParameters()
        {
            SpecFilePath = "";
            DatabaseFilePath = "";
            OutputDir = null;
            SearchModeInt = 1;
            InternalCleavageMode = InternalCleavageType.SingleInternalCleavage;
            TagBasedSearch = true;
            TdaInt = 0;
            TargetDecoySearchMode = DatabaseSearchMode.Target;
            PrecursorIonTolerancePpm = 10;
            ProductIonTolerancePpm = 10;
            MinSequenceLength = 21;
            MaxSequenceLength = 500;
            MinPrecursorIonCharge = 2;
            MaxPrecursorIonCharge = 50;
            MinProductIonCharge = 1;
            MaxProductIonCharge = 20;
            MinSequenceMass = 3000;
            MaxSequenceMass = 50000;
            FeatureFilePath = "";
            MaxNumThreads = 0;
            ActivationMethod = ActivationMethod.Unknown;
            ScansFilePath = "";

            // Single-source the defaults...
            var o = new MsPfParameters();
            InternalCleavageMode = o.InternalCleavageMode;
            TargetDecoySearchMode = o.TargetDecoySearchMode;
            TagBasedSearch = o.TagBasedSearch;
            PrecursorIonTolerancePpm = o.PrecursorIonTolerancePpm;
            ProductIonTolerancePpm = o.ProductIonTolerancePpm;
            MinSequenceLength = o.MinSequenceLength;
            MaxSequenceLength = o.MaxSequenceLength;
            MinPrecursorIonCharge = o.MinPrecursorIonCharge;
            MaxPrecursorIonCharge = o.MaxPrecursorIonCharge;
            MinProductIonCharge = o.MinProductIonCharge;
            MaxProductIonCharge = o.MaxProductIonCharge;
            MinSequenceMass = o.MinSequenceMass;
            MaxSequenceMass = o.MaxSequenceMass;
            ActivationMethod = o.ActivationMethod;
            MaxNumThreads = o.MaxNumThreads;
            // TODO: ...
        }

        /// <summary>
        /// Specific MS2 scans to process
        /// </summary>
        public SortedSet<int> ScanNumbers { get; set; }

        private IEnumerable<SearchModification> _searchModifications;
        private int _maxNumDynModsPerSequence;



        public bool Validate()
        {
            // Spec file path validation
            if (string.IsNullOrWhiteSpace(SpecFilePath))
            {
                PrintError("Missing parameter for spectrum file path");
                return false;
            }

            // Check for folder-type datasets, and replace specFilePath with the directory name if it is.
            SpecFilePath = MassSpecDataReaderFactory.GetDatasetName(SpecFilePath);

            var isDirectoryDataset = MassSpecDataReaderFactory.IsADirectoryDataset(SpecFilePath);
            // True if specFilePath is a directory that is NOT a supported folder-type dataset.
            var specPathIsDirectory = Directory.Exists(SpecFilePath) && !isDirectoryDataset;

            if (!File.Exists(SpecFilePath) && !specPathIsDirectory && !isDirectoryDataset)
            {
                PrintError("File not found: " + SpecFilePath);
                return false;
            }

            var types = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;

            if (!specPathIsDirectory && !(types.Select(ext => SpecFilePath.ToLower().EndsWith(ext)).Any()))
            {
                PrintError("Invalid file extension for spectrum file: (" + Path.GetExtension(SpecFilePath) + ") " + SpecFilePath);
                return false;
            }

            // TODO: Handle non-.raw files in the subfolder
            SpecFilePaths = Directory.Exists(SpecFilePath) && !MassSpecDataReaderFactory.IsADirectoryDataset(SpecFilePath) ? Directory.GetFiles(SpecFilePath, "*.raw") : new[] { SpecFilePath };

            // Database path validation
            if (string.IsNullOrWhiteSpace(DatabaseFilePath))
            {
                PrintError("Missing parameter for database file path");
                return false;
            }
            if (!File.Exists(DatabaseFilePath))
            {
                PrintError("File not found: " + DatabaseFilePath);
                return false;
            }
            var dbExtension = Path.GetExtension(DatabaseFilePath).ToLower();
            if (!dbExtension.Equals(".fa") &&
                !dbExtension.Equals(".fasta"))
            {
                PrintError("Invalid extension for the database file path (" + dbExtension + ")");
                return false;
            }

            // Output directory validation
            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                // Must use "Path.GetFullPath" to return the absolute path when the source file is in the working directory
                // But, it could cause problems with too-long paths.
                OutputDir = specPathIsDirectory ? SpecFilePath : Path.GetDirectoryName(Path.GetFullPath(SpecFilePath));
            }

            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                PrintError("Invalid output file directory: " + OutputDir);
                return false;
            }

            if (!Directory.Exists(OutputDir))
            {
                if (File.Exists(OutputDir) && !File.GetAttributes(OutputDir).HasFlag(FileAttributes.Directory))
                {
                    PrintError("OutputDir \"" + OutputDir + "\" is not a directory!");
                    return false;
                }
                Directory.CreateDirectory(OutputDir);
            }

            // Mods file validation
            if (!string.IsNullOrWhiteSpace(ModsFilePath) && !File.Exists(ModsFilePath))
            {
                PrintError("Modifications file not found: " + ModsFilePath);
                return false;
            }

            try
            {
                var errorMessage = LoadModsFile(ModsFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                {
                    PrintError(errorMessage);
                    return false;
                }
            }
            catch (Exception ex)
            {
                PrintError("Exception parsing the file for parameter -mod: " + ex.Message);
                return false;
            }

            // Scans file validation
            if (!string.IsNullOrWhiteSpace(ScansFilePath) && !File.Exists(ScansFilePath))
            {
                PrintError("Scans File file not found: " + ScansFilePath);
                return false;
            }
            try
            {
                var errorMessage = LoadScansFile(ScansFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                {
                    PrintError(errorMessage);
                    return false;
                }
            }
            catch (Exception ex)
            {
                PrintError("Exception parsing the file for parameter -scansFile: " + ex.Message);
                return false;
            }

            // Feature file validation
            if (!string.IsNullOrWhiteSpace(FeatureFilePath) && !File.Exists(FeatureFilePath))
            {
                PrintError("Feature File not found: " + FeatureFilePath);
                return false;
            }
            if (!string.IsNullOrWhiteSpace(FeatureFilePath) &&
                !Path.GetExtension(FeatureFilePath).ToLower().Equals(".csv") &&
                !Path.GetExtension(FeatureFilePath).ToLower().Equals(".ms1ft") &&
                !Path.GetExtension(FeatureFilePath).ToLower().Equals(".msalign"))
            {
                PrintError("Invalid extension for the Feature file path (" + Path.GetExtension(FeatureFilePath) + ")");
                return false;
            }

            // MinX/MaxX validation
            if (MinSequenceLength > MaxSequenceLength)
            {
                PrintError("MinPrecursorCharge (" + MinPrecursorIonCharge + ") is larger than MaxPrecursorCharge (" + MaxPrecursorIonCharge + ")!");
                return false;
            }

            if (MinProductIonCharge > MaxProductIonCharge)
            {
                PrintError("MinFragmentCharge (" + MinProductIonCharge + ") is larger than MaxFragmentCharge (" + MaxProductIonCharge + ")!");
                return false;
            }

            if (MinSequenceMass > MaxSequenceMass)
            {
                PrintError("MinSequenceMassInDa (" + MinSequenceMass + ") is larger than MaxSequenceMassInDa (" + MaxSequenceMass + ")!");
                return false;
            }

            MaxNumThreads = GetOptimalMaxThreads(MaxNumThreads);

            return true;
        }

        private static void PrintError(string errorMessage)
        {
            Console.WriteLine();
            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine("Error: " + errorMessage);
            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine();
        }

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
            Console.WriteLine("InternalCleavageMode: " + InternalCleavageMode);
            Console.WriteLine("Tag-based search: " + TagBasedSearch);
            Console.WriteLine("Tda: " + (TargetDecoySearchMode == DatabaseSearchMode.Both ? "Target+Decoy" : TargetDecoySearchMode.ToString()));
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
                    writer.WriteLine("InternalCleavageMode\t" + InternalCleavageMode);
                    writer.WriteLine("Tag-based search\t" + TagBasedSearch);
                    writer.WriteLine("Tda\t" + (TargetDecoySearchMode == DatabaseSearchMode.Both ? "Target+Decoy" : TargetDecoySearchMode.ToString()));
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
