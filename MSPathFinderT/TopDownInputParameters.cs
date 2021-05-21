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
using PRISM.AppSettings;
using PRISM.Logging;

namespace MSPathFinderT
{
    public class TopDownInputParameters : MsPfParameters
    {
        // Ignore Spelling: Da, frag, hyperthreaded, ic, tda, tol, wildcards

        [Option("i", "s", "specFile", ArgPosition = 1, Required = true,
            HelpText = "Spectrum File (.raw or .pbf)",
            HelpShowsDefault = false)]
        public override string SpecFilePath { get; set; }

        public List<FileSystemInfo> SpecFilePaths { get; }

        [Option("d", "database", Required = true,
            HelpText = "Database File (*.fasta or *.fa or *.faa)",
            HelpShowsDefault = false)]
        public override string DatabaseFilePath { get; set; }

        [Option("o", "outputDir",
            HelpText = "Output Directory",
            HelpShowsDefault = false)]
        public override string OutputDir { get; set; }

        [Option("m", "searchMode", Min = 0, Max = 2,
            HelpText = "Search Mode (old format) (0: multiple internal cleavages, 1: single internal cleavage, 2: no internal cleavage)",
            Hidden = true)]
        [Obsolete("Use InternalCleavageMode")]
        public override int SearchModeInt
        {
            get
            {
                if (InternalCleavageMode == InternalCleavageType.MultipleInternalCleavages)
                {
                    return 0;
                }

                if (InternalCleavageMode == InternalCleavageType.SingleInternalCleavage)
                {
                    return 1;
                }

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

        [Option("ic",
            HelpText = "Search Mode")]
        public override InternalCleavageType InternalCleavageMode { get; set; }

        [Option("tagSearch",
            HelpText = "Include Tag-based Search (use true or false;\nor use '0' for false or '1' for true)")]
        public override bool TagBasedSearch { get; set; }

        [Option("memMatches", "MatchesPerSpectrumToKeepInMemory",
            HelpText = "Number of matches to keep in memory; these matches are used when computing spectral E-values")]
        public override int MatchesPerSpectrumToKeepInMemory { get; set; }

        [Option("n", "NumMatchesPerSpec", "MatchesPerSpectrumToReport",
            HelpText = "Number of results to report for each mass spectrum")]
        public override int MatchesPerSpectrumToReport { get; set; }

        [Option("IncludeDecoy", "IncludeDecoys", "IncludeDecoyResults",
            HelpText = "Include decoy results in the _IcTda.tsv file")]
        public override bool IncludeDecoyResults { get; set; }

        [Option("mod",
            HelpText = "Path to modification file that defines static and dynamic modifications. " +
                       "Modifications can alternatively be defined in a parameter file, as specified by /ParamFile or -ParamFile\n" +
                       "Modifications defined using the -mod switch take precedence over modifications defined in a parameter file\n" +
                       "(Default: empty string, meaning no modifications)",
            HelpShowsDefault = false
            )]
        public string ModsFilePath { get; set; }

        // ReSharper disable once UnusedMember.Global
        [Option("tda", Min = -1, Max = 1,
            HelpText = "Database search mode:\n" +
                       "0: don't search decoy database, \n" +
                       "1: search shuffled decoy database\n")]
        public int TdaInt
        {
            get
            {
                switch (TargetDecoySearchMode)
                {
                    case DatabaseSearchMode.Both:
                        // Target and Decoy search
                        return 1;

                    case DatabaseSearchMode.Decoy:
                        // Decoy search only (shuffled database)
                        return -1;

                    default:
                        // Target search only
                        return 0;
                }
            }
            set
            {
                switch (value)
                {
                    case -1:
                        TargetDecoySearchMode = DatabaseSearchMode.Decoy;
                        break;
                    case 1:
                        TargetDecoySearchMode = DatabaseSearchMode.Both;
                        break;
                    default:
                        TargetDecoySearchMode = DatabaseSearchMode.Target;
                        break;
                }
            }
        }

        /// <summary>
        /// Database search mode enum
        /// </summary>
        /// <remarks>This can alternatively be set using TdaInt</remarks>
        public override DatabaseSearchMode TargetDecoySearchMode { get; set; }

        [Option("overwrite", "OverwriteExistingResults",
            HelpText = "Overwrite existing results. If false (default), looks for files _IcTarget.tsv and _IcDecoy.tsv and uses the results in the files if found")]
        public override bool OverwriteExistingResults { get; set; }

        [Option("t", "precursorTol", "PMTolerance", /*Min = 1,*/
            HelpText = "Precursor Tolerance (in PPM)")]
        public override double PrecursorIonTolerancePpm
        {
            get => PrecursorIonTolerance.GetValue();
            set => PrecursorIonTolerance = new Tolerance(value);
        }

        [Option("f", "fragmentTol", "FragTolerance", /*Min = 1,*/
            HelpText = "Fragment Ion Tolerance (in PPM)")]
        public override double ProductIonTolerancePpm
        {
            get => ProductIonTolerance.GetValue();
            set => ProductIonTolerance = new Tolerance(value);
        }

        [Option("minLength", Min = 0,
            HelpText = "Minimum Sequence Length")]
        public override int MinSequenceLength { get; set; }

        [Option("maxLength", Min = 0,
            HelpText = "Maximum Sequence Length")]
        public override int MaxSequenceLength { get; set; }

        [Option("minCharge", Min = 1,
            HelpText = "Minimum precursor ion charge")]
        public override int MinPrecursorIonCharge { get; set; }

        [Option("maxCharge", Min = 1,
            HelpText = "Maximum precursor ion charge")]
        public override int MaxPrecursorIonCharge { get; set; }

        [Option("minFragCharge", Min = 1,
            HelpText = "Minimum fragment ion charge")]
        public override int MinProductIonCharge { get; set; }

        [Option("maxFragCharge", Min = 1,
            HelpText = "Maximum fragment ion charge")]
        public override int MaxProductIonCharge { get; set; }

        [Option("minMass", /*Min = 1,*/
            HelpText = "Minimum sequence mass in Da")]
        public override double MinSequenceMass { get; set; }

        [Option("maxMass", /*Min = 1,*/
            HelpText = "Maximum sequence mass in Da")]
        public override double MaxSequenceMass { get; set; }

        [Option("feature",
            HelpText = ".ms1ft, _isos.csv, or .msalign feature file (typically the results from ProMex); " +
                       "leave blank/undefined if processing multiple input files",
            HelpShowsDefault = false)]
        public override string FeatureFilePath { get; set; }

        [Option("threads", Min = 0,
            HelpText = "Maximum number of threads, or 0 to set automatically")]
        public override int MaxNumThreads { get; set; }

        [Option("act", "ActivationMethod",
            HelpText = "Activation Method")]
        public override ActivationMethod ActivationMethod { get; set; }

        [Option("scansFile", "ScansFilePath",
            HelpText = "Optional text file with MS2 scans to process (tab, comma, or space separated); " +
                       "any integer in the file is assumed to be a scan number to process",
            HelpShowsDefault = false)]
        public string ScansFilePath { get; set; }

        [Option("flip",
            HelpText = "If specified, FLIP scoring code will be used\n(supports UVPD spectra)")]
        public bool UseFLIP { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The CommandLineParser class will look for a parameter named /ParamFile (or -ParamFile)
        /// If defined, it will read settings from that parameter file, looking for settings that match the properties in this class
        /// It will not read NumMods, StaticMod, or DynamicMod entries
        /// Reading that info is accomplished via methods in this class
        /// </remarks>
        public TopDownInputParameters()
        {
            ScansFilePath = string.Empty;
            UseFLIP = false;

            SpecFilePaths = new List<FileSystemInfo>();
        }

        /// <summary>
        /// Validate options
        /// </summary>
        /// <parameterFilePath>Path to the parameter file specified by /ParamFile or -ParamFile, or an empty string if not provided</parameterFilePath>
        /// <returns>True if options are valid, otherwise false</returns>
        /// <remarks>Will also load dynamic and static modifications from the parameter file, if defined</remarks>
        public bool Validate(string parameterFilePath)
        {
            var defaultOutputDirectoryPath = string.Empty;
            SpecFilePaths.Clear();

            // Spec file path validation
            if (string.IsNullOrWhiteSpace(SpecFilePath))
            {
                ShowError("Missing parameter for spectrum file path");
                return false;
            }

            if (SpecFilePath.Contains("*") || SpecFilePath.Contains("?"))
            {
                // SpecFilePath has wildcards
                // Validate each matching file or directory

                var cleanPath = SpecFilePath.Replace('*', '_').Replace('?', '_');
                var lastSlash = SpecFilePath.LastIndexOf(Path.DirectorySeparatorChar);

                string wildcardSpec;
                if (lastSlash >= 0)
                {
                    wildcardSpec = SpecFilePath.Substring(lastSlash + 1);
                }
                else
                {
                    wildcardSpec = SpecFilePath;
                }

                var sourcePathInfo = new FileInfo(cleanPath);
                var workingDirectory = sourcePathInfo.Directory ?? new DirectoryInfo(".");

                var matchingFiles = workingDirectory.GetFiles(wildcardSpec).ToList();
                if (matchingFiles.Count > 0)
                {
                    Console.WriteLine("Finding dataset files that match " + SpecFilePath);
                    foreach (var datasetFile in matchingFiles)
                    {
                        ShowDebug(PathUtils.CompactPathString(datasetFile.FullName, 80), 0);
                        var specFiles = ValidateSourceData(datasetFile.FullName, out _);
                        if (specFiles.Count == 0)
                            return false;

                        SpecFilePaths.AddRange(specFiles);
                        defaultOutputDirectoryPath = workingDirectory.FullName;
                    }

                    Console.WriteLine();
                }
                else
                {
                    var matchingDirectories = workingDirectory.GetDirectories(wildcardSpec).ToList();
                    if (matchingDirectories.Count > 0)
                    {
                        Console.WriteLine("Finding dataset directories that match " + SpecFilePath);
                        foreach (var datasetDirectory in matchingDirectories)
                        {
                            ShowDebug(PathUtils.CompactPathString(datasetDirectory.FullName, 80), 0);
                            var specFiles = ValidateSourceData(datasetDirectory.FullName, out _);
                            if (specFiles.Count == 0)
                                return false;

                            SpecFilePaths.AddRange(specFiles);
                            defaultOutputDirectoryPath = workingDirectory.FullName;
                        }

                        Console.WriteLine();
                    }
                    else
                    {
                        ShowWarning(string.Format(
                            "No files or directories matched '{0}' in directory {1}", wildcardSpec, workingDirectory.FullName));

                        return false;
                    }
                }
            }
            else
            {
                var specFiles = ValidateSourceData(SpecFilePath, out var specPathIsDirectory);
                if (specFiles.Count == 0)
                    return false;

                SpecFilePaths.AddRange(specFiles);

                if (specPathIsDirectory)
                {
                    var datasetDirectory = new DirectoryInfo(SpecFilePath);
                    defaultOutputDirectoryPath = datasetDirectory.FullName;
                }
                else
                {
                    var datasetFile = new FileInfo(SpecFilePath);
                    defaultOutputDirectoryPath = datasetFile.DirectoryName ?? ".";
                }
            }

            // Database path validation
            if (string.IsNullOrWhiteSpace(DatabaseFilePath))
            {
                ShowError("Missing parameter for database file path");
                return false;
            }
            if (!File.Exists(DatabaseFilePath))
            {
                ShowError("File not found: " + DatabaseFilePath);
                return false;
            }

            if (!FastaDatabaseConstants.ValidFASTAExtension(DatabaseFilePath))
            {
                ShowError("Invalid extension for the database file path (" + Path.GetExtension(DatabaseFilePath) + ")");
                return false;
            }

            // Output directory validation
            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                OutputDir = defaultOutputDirectoryPath;
            }

            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                ShowError("Invalid output file directory: " + OutputDir);
                return false;
            }

            if (!Directory.Exists(OutputDir))
            {
                if (File.Exists(OutputDir) && !File.GetAttributes(OutputDir).HasFlag(FileAttributes.Directory))
                {
                    ShowError("OutputDir \"" + OutputDir + "\" is not a directory!");
                    return false;
                }
                Directory.CreateDirectory(OutputDir);
            }

            // Mods file validation
            if (string.IsNullOrWhiteSpace(ModsFilePath))
            {
                if (string.IsNullOrWhiteSpace(parameterFilePath))
                {
                    AminoAcidSet = new AminoAcidSet();
                    Modifications = new List<SearchModification>();
                }
                else
                {
                    var success = LoadModsFromMSPathFinderParameterFile(parameterFilePath);
                    if (!success)
                    {
                        return false;
                    }
                }
            }
            else
            {
                if (!File.Exists(ModsFilePath))
                {
                    ShowError("Modifications file not found: " + ModsFilePath);
                    return false;
                }

                var errorMessage = LoadModsFile(ModsFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                {
                    ShowError(errorMessage);
                    return false;
                }
            }

            // Scans file validation
            if (!string.IsNullOrWhiteSpace(ScansFilePath) && !File.Exists(ScansFilePath))
            {
                ShowError("Scans File not found: " + ScansFilePath);
                return false;
            }
            try
            {
                var errorMessage = LoadScansFile(ScansFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                {
                    ShowError(errorMessage);
                    return false;
                }
            }
            catch (Exception ex)
            {
                ShowError("Exception parsing the file for parameter -scansFile: " + ex.Message);
                return false;
            }

            // Feature file validation
            if (!string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                if (FeatureFilePath.Contains("*") || FeatureFilePath.Contains("?"))
                {
                    // Feature file path has a wildcard; change to an empty string so that the file is auto-found
                    FeatureFilePath = string.Empty;
                }
                else if (!File.Exists(FeatureFilePath))
                {
                    ShowError("Feature File not found: " + FeatureFilePath);
                    return false;
                }
                else if (SpecFilePaths.Count > 1)
                {
                    ShowWarning(string.Format(
                        "Processing multiple spectrum files, but a single feature file is defined:\n  {0}\n\n" +
                        "You probably want to leave the feature file path undefined\n" +
                        "by either not using the -feature switch or\n" +
                        "by not including 'feature=FilePath' in the parameter file",
                        FeatureFilePath));
                }
            }

            if (!string.IsNullOrWhiteSpace(FeatureFilePath) &&
                !Path.GetExtension(FeatureFilePath).Equals(".csv", StringComparison.OrdinalIgnoreCase) &&
                !Path.GetExtension(FeatureFilePath).Equals(".ms1ft", StringComparison.OrdinalIgnoreCase) &&
                !Path.GetExtension(FeatureFilePath).Equals(".msalign", StringComparison.OrdinalIgnoreCase))
            {
                ShowError("Invalid extension for the Feature file path (" + Path.GetExtension(FeatureFilePath) + ")");
                return false;
            }

            if (MinSequenceLength > MaxSequenceLength)
            {
                ShowError("MinPrecursorCharge (" + MinPrecursorIonCharge + ") is larger than MaxPrecursorCharge (" + MaxPrecursorIonCharge + ")!");
                return false;
            }

            if (MinProductIonCharge > MaxProductIonCharge)
            {
                ShowError("MinFragmentCharge (" + MinProductIonCharge + ") is larger than MaxFragmentCharge (" + MaxProductIonCharge + ")!");
                return false;
            }

            if (MinSequenceMass > MaxSequenceMass)
            {
                ShowError("MinSequenceMassInDa (" + MinSequenceMass + ") is larger than MaxSequenceMassInDa (" + MaxSequenceMass + ")!");
                return false;
            }

            MaxNumThreads = GetOptimalMaxThreads(MaxNumThreads);

            if (MatchesPerSpectrumToReport >= MatchesPerSpectrumToKeepInMemory)
            {
                MatchesPerSpectrumToKeepInMemory = MatchesPerSpectrumToReport + 1;
            }

            return true;
        }

        /// <summary>
        /// Validate that the dataset is a valid mass spec dataset
        /// </summary>
        /// <param name="sourceDatasetPath">Path to a dataset file or directory</param>
        /// <param name="specPathIsDirectory">Output: true if the specified path is a directory</param>
        /// <returns>List of files to process; empty list if an error</returns>
        private List<FileSystemInfo> ValidateSourceData(string sourceDatasetPath, out bool specPathIsDirectory)
        {
            // Check for folder-type datasets, and replace sourceDatasetPath with the directory name if it is.
            sourceDatasetPath = MassSpecDataReaderFactory.GetDatasetName(sourceDatasetPath);

            var isDirectoryDataset = MassSpecDataReaderFactory.IsADirectoryDataset(sourceDatasetPath);

            // This will be True if sourceDatasetPath is a directory that is NOT a supported folder-type dataset.
            specPathIsDirectory = Directory.Exists(sourceDatasetPath) && !isDirectoryDataset;

            if (!File.Exists(sourceDatasetPath) && !specPathIsDirectory && !isDirectoryDataset)
            {
                ShowError("File not found: " + sourceDatasetPath);
                return new List<FileSystemInfo>();
            }

            // The extensions in this variable start with a period
            var supportedFileExtensions = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;

            if (!specPathIsDirectory && !supportedFileExtensions.Select(ext => sourceDatasetPath.ToLower().EndsWith(ext)).Any())
            {
                ShowError("Invalid file extension for spectrum file (" + Path.GetExtension(sourceDatasetPath) + "): " + sourceDatasetPath);
                return new List<FileSystemInfo>();
            }

            var specFiles = GetSpecFilePaths(sourceDatasetPath).ToList();
            return specFiles;
        }

        private static void ShowDebug(string message, int emptyLinesBeforeMessage = 1)
        {
            ConsoleMsgUtils.ShowDebugCustom(message, emptyLinesBeforeMessage: emptyLinesBeforeMessage);
        }

        private static void ShowError(string errorMessage, Exception ex = null)
        {
            ShowErrorOrWarning(errorMessage, ex);
        }

        private static void ShowWarning(string message)
        {
            ShowErrorOrWarning(message, null, string.Empty);
        }

        private static void ShowErrorOrWarning(string message, Exception ex = null, string messagePrefix = "Error: ")
        {
            Console.WriteLine();
            ConsoleMsgUtils.ShowWarning("{0}\n{1}{2}\n{0}",
                "----------------------------------------------------------",
                messagePrefix, message);

            if (ex != null)
            {
                ConsoleMsgUtils.ShowWarning(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
            }

            Console.WriteLine();
        }

        public void Display(string parameterFilePath)
        {
            foreach (var specFilePath in SpecFilePaths)
            {
                Console.WriteLine("{0,-28} {1}", "SpectrumFilePath:", specFilePath);
            }

            Console.WriteLine("{0,-28} {1}", "DatabaseFilePath:", DatabaseFilePath);

            var featureFilePathToShow = string.IsNullOrWhiteSpace(FeatureFilePath) ? "auto-find" : FeatureFilePath;
            Console.WriteLine("{0,-28} {1}", "FeatureFilePath:", featureFilePathToShow);

            var parameterFilePathToShow = string.IsNullOrWhiteSpace(parameterFilePath) ? "N/A" : parameterFilePath;
            Console.WriteLine("{0,-28} {1}", "ParameterFilePath:", parameterFilePathToShow);
            Console.WriteLine("{0,-28} {1}", "OutputDir:", OutputDir);
            Console.WriteLine();

            Console.WriteLine("{0,-29} {1}", "MaxThreads:", MaxNumThreads);
            Console.WriteLine("{0,-29} {1}", "InternalCleavageMode:", InternalCleavageMode);
            Console.WriteLine("{0,-29} {1}", "Tag-based search:", TagBasedSearch);
            Console.WriteLine("{0,-29} {1}", "Tda:", TargetDecoySearchMode == DatabaseSearchMode.Both ? "Target+Decoy" : TargetDecoySearchMode.ToString());
            Console.WriteLine("{0,-29} {1}", "PrecursorIonTolerancePpm:", PrecursorIonTolerancePpm);
            Console.WriteLine("{0,-29} {1}", "ProductIonTolerancePpm:", ProductIonTolerancePpm);
            Console.WriteLine("{0,-29} {1}", "MinSequenceLength:", MinSequenceLength);
            Console.WriteLine("{0,-29} {1}", "MaxSequenceLength:", MaxSequenceLength);
            Console.WriteLine("{0,-29} {1}", "MinPrecursorIonCharge:", MinPrecursorIonCharge);
            Console.WriteLine("{0,-29} {1}", "MaxPrecursorIonCharge:", MaxPrecursorIonCharge);
            Console.WriteLine("{0,-29} {1}", "MinProductIonCharge:", MinProductIonCharge);
            Console.WriteLine("{0,-29} {1}", "MaxProductIonCharge:", MaxProductIonCharge);
            Console.WriteLine("{0,-29} {1}", "MinSequenceMass:", MinSequenceMass);
            Console.WriteLine("{0,-29} {1}", "MaxSequenceMass:", MaxSequenceMass);
            Console.WriteLine("{0,-29} {1}", "MatchesPerSpectrumToReport:", MatchesPerSpectrumToReport);
            Console.WriteLine("{0,-29} {1}", "MatchesPerSpecToKeepInMemory:", MatchesPerSpectrumToKeepInMemory);
            Console.WriteLine("{0,-29} {1}", "IncludeDecoyResults:", IncludeDecoyResults);
            Console.WriteLine("{0,-29} {1}", "MaxDynamicModsPerSequence:", MaxDynamicModificationsPerSequence);
            Console.WriteLine("{0,-29} {1}", "OverwriteExistingResults:", OverwriteExistingResults);

            Console.WriteLine();
            Console.WriteLine("Modifications:");

            foreach (var searchMod in Modifications)
            {
                Console.WriteLine(searchMod.ToString(true));
            }

            if (!string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                Console.WriteLine();
                Console.WriteLine("Loading MS1 features from " + FeatureFilePath);
            }

            if (ScanNumbers?.Any() == true)
            {
                Console.WriteLine("Processing specific MS2 scans:");
                Console.WriteLine(string.Join(", ", ScanNumbers));
            }

            if (UseFLIP)
            {
                Console.WriteLine("Using FLIP scoring.");
            }

            Console.WriteLine();
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

        private static IEnumerable<FileSystemInfo> GetSpecFilePaths(string fileOrDirectoryPath)
        {
            try
            {
                if (!Directory.Exists(fileOrDirectoryPath))
                {
                    // Assume a file
                    return new[] { new FileInfo(fileOrDirectoryPath) };
                }

                var workingDirectory = new DirectoryInfo(fileOrDirectoryPath);

                if (MassSpecDataReaderFactory.IsADirectoryDataset(fileOrDirectoryPath))
                {
                    return new[] { workingDirectory };
                }

                // Process all .raw files in the given directory
                return workingDirectory.GetFiles("*.raw");
            }
            catch (Exception ex)
            {
                ShowError(string.Format("Exception examining the file or directory path ({0}): {1}", fileOrDirectoryPath, ex.Message), ex);
                return new List<FileInfo>();
            }
        }

        /// <summary>
        /// Parse the modification info file, if defined
        /// </summary>
        /// <param name="modFilePath"></param>
        /// <returns>Error message if an error</returns>
        private string LoadModsFile(string modFilePath)
        {
            try
            {
                if (!File.Exists(modFilePath))
                {
                    return "-mod file not found: " + modFilePath;
                }

                var parser = new ModFileParser(modFilePath);
                Modifications = parser.SearchModifications.ToList();
                MaxDynamicModificationsPerSequence = parser.MaxNumDynModsPerSequence;

                if (Modifications == null)
                {
                    return "Error while parsing " + modFilePath;
                }

                AminoAcidSet = new AminoAcidSet(Modifications, MaxDynamicModificationsPerSequence);
                return string.Empty;
            }
            catch (Exception ex)
            {
                return "Exception parsing the file for parameter -mod: " + ex.Message;
            }
        }

        private string LoadScansFile(string scansFilePath)
        {
            if (string.IsNullOrWhiteSpace(scansFilePath))
            {
                return string.Empty;
            }

            if (!File.Exists(scansFilePath))
            {
                return "-scansFile not found: " + scansFilePath;
            }

            var scanNumbers = new SortedSet<int>();

            var delimiters = new[] { ' ', '\t', ',' };

            using (var reader = new StreamReader(new FileStream(scansFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                while (!reader.EndOfStream)
                {
                    var dataLine = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(dataLine))
                    {
                        continue;
                    }

                    var dataValues = dataLine.Split(delimiters);
                    foreach (var value in dataValues)
                    {
                        if (!int.TryParse(value, out var scanNumber))
                        {
                            continue;
                        }

                        if (!scanNumbers.Contains(scanNumber))
                        {
                            scanNumbers.Add(scanNumber);
                        }
                    }
                }
            }

            ScanNumbers = scanNumbers;

            return string.Empty;
        }

        /// <summary>
        /// Read the MSPathFinder parameter file to look for NumMods plus any StaticMod or DynamicMod entries
        /// </summary>
        /// <param name="parameterFilePath"></param>
        /// <returns>True if success, false if an error</returns>
        private bool LoadModsFromMSPathFinderParameterFile(string parameterFilePath)
        {
            var paramFileReader = new KeyValueParamFileReader("MSPathFinder", parameterFilePath);
            RegisterEvents(paramFileReader);

            var success = paramFileReader.ParseKeyValueParameterFile(out var paramFileEntries);
            if (!success)
            {
                ConsoleMsgUtils.ShowWarning(paramFileReader.ErrorMessage);
                return false;
            }

            var numMods = 0;
            var staticMods = new List<string>();
            var dynamicMods = new List<string>();

            try
            {
                foreach (var kvSetting in paramFileEntries)
                {
                    var paramValue = kvSetting.Value;

                    if (kvSetting.Key.Equals("NumMods", StringComparison.OrdinalIgnoreCase))
                    {
                        if (int.TryParse(paramValue, out var intValue))
                        {
                            numMods = intValue;
                        }
                        else
                        {
                            const string errMsg = "Invalid value for NumMods in the MSPathFinder parameter file";
                            ConsoleMsgUtils.ShowWarning(errMsg + ": " + kvSetting.Key + "=" + kvSetting.Value);
                            return false;
                        }
                    }
                    else if (kvSetting.Key.Equals("StaticMod", StringComparison.OrdinalIgnoreCase))
                    {
                        if (!string.IsNullOrWhiteSpace(paramValue) && !paramValue.Equals("none", StringComparison.OrdinalIgnoreCase))
                        {
                            staticMods.Add(paramValue);
                        }
                    }
                    else if (kvSetting.Key.Equals("DynamicMod", StringComparison.OrdinalIgnoreCase))
                    {
                        if (!string.IsNullOrWhiteSpace(paramValue) && !paramValue.Equals("none", StringComparison.OrdinalIgnoreCase))
                        {
                            dynamicMods.Add(paramValue);
                        }
                    }
                }

                var modsValidated = StoreMSPathFinderModifications(numMods, staticMods, dynamicMods);

                return modsValidated;
            }
            catch (Exception ex)
            {
                ShowError("Exception extracting dynamic and static mod information from the MSPathFinder parameter file", ex);
                return false;
            }
        }

        private bool StoreMSPathFinderModifications(int maxDynamicModsPerSequence, IEnumerable<string> staticMods, IEnumerable<string> dynamicMods)
        {
            try
            {
                MaxDynamicModificationsPerSequence = maxDynamicModsPerSequence;

                Modifications = new List<SearchModification>();

                foreach (var staticMod in staticMods)
                {
                    if (ParseMSPathFinderValidateMod(staticMod, out var modClean, out var parsedMods))
                    {
                        if (modClean.Contains(",opt,"))
                        {
                            // Static (fixed) mod is listed as dynamic
                            // Abort the analysis since the parameter file is misleading and needs to be fixed
                            const string errMsg = "Static mod definition contains ',opt,'; update the param file to have ',fix,' or change to 'DynamicMod='";
                            ShowError(errMsg + "\n" + staticMod);
                            return false;
                        }

                        Modifications.AddRange(parsedMods);
                    }
                    else
                    {
                        return false;
                    }
                }

                foreach (var dynamicMod in dynamicMods)
                {
                    if (ParseMSPathFinderValidateMod(dynamicMod, out var modClean, out var parsedMods))
                    {
                        if (modClean.Contains(",fix,"))
                        {
                            // Dynamic (optional) mod is listed as static
                            // Abort the analysis since the parameter file is misleading and needs to be fixed
                            const string errMsg = "Dynamic mod definition contains ',fix,'; update the param file to have ',opt,' or change to 'StaticMod='";
                            ShowError(errMsg + "\n" + dynamicMod);
                            return false;
                        }

                        Modifications.AddRange(parsedMods);
                    }
                    else
                    {
                        return false;
                    }
                }

                AminoAcidSet = new AminoAcidSet(Modifications, MaxDynamicModificationsPerSequence);

                return true;
            }
            catch (Exception ex)
            {
                ShowError("Exception storing dynamic and static mod information loaded from the MSPathFinder parameter file", ex);
                return false;
            }
        }

        /// <summary>
        /// Validates that the modification definition text
        /// </summary>
        /// <param name="mod">Modification definition</param>
        /// <param name="modClean">Cleaned-up modification definition (output param)</param>
        /// <param name="parsedMods">
        /// List of search modifications determined using modClean
        /// There will be multiple items if the modification affects multiple residues
        /// </param>
        /// <returns>True if valid; false if invalid</returns>
        /// <remarks>A valid modification definition contains 5 parts and doesn't contain any whitespace</remarks>
        private bool ParseMSPathFinderValidateMod(string mod, out string modClean, out List<SearchModification> parsedMods)
        {
            var comment = string.Empty;

            modClean = string.Empty;

            var poundIndex = mod.IndexOf('#');
            if (poundIndex > 0)
            {
                comment = mod.Substring(poundIndex);
                mod = mod.Substring(0, poundIndex - 1).Trim();
            }

            var splitMod = mod.Split(',');

            if (splitMod.Length < 5)
            {
                // Invalid mod definition; must have 5 sections
                ConsoleMsgUtils.ShowWarning("Invalid modification string; must have 5 sections: " + mod);
                parsedMods = new List<SearchModification>();
                return false;
            }

            // Make sure mod does not have both * and any
            if (splitMod[1].Trim() == "*" && splitMod[3].ToLower().Trim() == "any")
            {
                ConsoleMsgUtils.ShowWarning("Modification cannot contain both * and any: " + mod);
                parsedMods = new List<SearchModification>();
                return false;
            }

            // Reconstruct the mod definition, making sure there is no whitespace
            modClean = splitMod[0].Trim();
            for (var index = 1; index <= splitMod.Length - 1; index++)
            {
                modClean += "," + splitMod[index].Trim();
            }

            if (!string.IsNullOrWhiteSpace(comment))
            {
                modClean += "     " + comment;
            }

            try
            {
                parsedMods = ModFileParser.ParseModification(modClean);
            }
            catch (Exception ex)
            {
                ConsoleMsgUtils.ShowWarning(
                    "Invalid modification line: {0}\n{1}\n\n(originally {2})",
                    modClean, ex.Message, mod);
                parsedMods = new List<SearchModification>();
                return false;
            }

            return parsedMods != null;
        }

        #region "Events"

        private static void RegisterEvents(IEventNotifier processingClass, bool writeDebugEventsToLog = true)
        {
            if (writeDebugEventsToLog)
            {
                processingClass.DebugEvent += DebugEventHandler;
            }
            else
            {
                processingClass.DebugEvent += DebugEventHandlerConsoleOnly;
            }

            processingClass.StatusEvent += StatusEventHandler;
            processingClass.ErrorEvent += ErrorEventHandler;
            processingClass.WarningEvent += WarningEventHandler;
            // Ignore: processingClass.ProgressUpdate += ProgressUpdateHandler;
        }

        private static void DebugEventHandlerConsoleOnly(string statusMessage)
        {
            ConsoleMsgUtils.ShowDebug(statusMessage);
        }

        private static void DebugEventHandler(string statusMessage)
        {
            ConsoleMsgUtils.ShowDebug(statusMessage);
        }

        private static void StatusEventHandler(string statusMessage)
        {
            Console.WriteLine(statusMessage);
        }

        private static void ErrorEventHandler(string errorMessage, Exception ex)
        {
            ShowError(errorMessage, ex);
        }

        private static void WarningEventHandler(string warningMessage)
        {
            ConsoleMsgUtils.ShowWarning(warningMessage);
        }

        #endregion
    }
}
