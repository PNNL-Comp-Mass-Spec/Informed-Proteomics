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
        // Ignore Spelling: ic, tda, tol, frag, Da

        [Option("s", "specFile", ArgPosition = 1, Required = true, HelpText = "Spectrum File (*.raw)", HelpShowsDefault = false)]
        public override string SpecFilePath { get; set; }

        public IEnumerable<string> SpecFilePaths { get; set; }

        [Option("d", "database", Required = true, HelpText = "Database File (*.fasta or *.fa or *.faa)", HelpShowsDefault = false)]
        public override string DatabaseFilePath { get; set; }

        [Option("o", "outputDir", HelpText = "Output Directory", HelpShowsDefault = false)]
        public override string OutputDir { get; set; }

        [Option("m", "searchMode", Min = 0, Max = 2, HelpText = "Search Mode (old format) (0: multiple internal cleavages, 1: single internal cleavage, 2: no internal cleavage)", Hidden = true)]
        [Obsolete("Use InternalCleavageMode")]
        public override int SearchModeInt
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
        public override InternalCleavageType InternalCleavageMode { get; set; }

        [Option("tagSearch", HelpText = "Include Tag-based Search (use true or false;\nor use '0' for false or '1' for true)")]
        public override bool TagBasedSearch { get; set; }

        [Option("n", "NumMatchesPerSpec", "MatchesPerSpectrumToReport", HelpText = "Number of results to report for each mass spectrum")]
        public override int MatchesPerSpectrumToReport { get; set; }

        [Option("IncludeDecoy", "IncludeDecoys", "IncludeDecoyResults", HelpText = "Include decoy results in the _IcTda.tsv file")]
        public override bool IncludeDecoyResults { get; set; }

        [Option("mod",
            HelpText = "Path to modification file that defines static and dynamic modifications. " +
                       "Modifications can alternatively be defined in a parameter file, as specified by /ParamFile or -ParamFile\n" +
                       "Modifications defined using the -mod switch take precedence over modifications defined in a parameter file\n" +
                       "(Default: empty string, meaning no modifications)",
            HelpShowsDefault = false)]
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

        [Option("t", "precursorTol", "PMTolerance", /*Min = 1,*/ HelpText = "Precursor Tolerance (in PPM)")]
        public override double PrecursorIonTolerancePpm
        {
            get => PrecursorIonTolerance.GetValue();
            set => PrecursorIonTolerance = new Tolerance(value);
        }

        [Option("f", "fragmentTol", "FragTolerance", /*Min = 1,*/ HelpText = "Fragment Ion Tolerance (in PPM)")]
        public override double ProductIonTolerancePpm
        {
            get => ProductIonTolerance.GetValue();
            set => ProductIonTolerance = new Tolerance(value);
        }

        [Option("minLength", Min = 0, HelpText = "Minimum Sequence Length")]
        public override int MinSequenceLength { get; set; }

        [Option("maxLength", Min = 0, HelpText = "Maximum Sequence Length")]
        public override int MaxSequenceLength { get; set; }

        [Option("minCharge", Min = 1, HelpText = "Minimum precursor ion charge")]
        public override int MinPrecursorIonCharge { get; set; }

        [Option("maxCharge", Min = 1, HelpText = "Maximum precursor ion charge")]
        public override int MaxPrecursorIonCharge { get; set; }

        [Option("minFragCharge", Min = 1, HelpText = "Minimum fragment ion charge")]
        public override int MinProductIonCharge { get; set; }

        [Option("maxFragCharge", Min = 1, HelpText = "Maximum fragment ion charge")]
        public override int MaxProductIonCharge { get; set; }

        [Option("minMass", /*Min = 1,*/ HelpText = "Minimum sequence mass in Da")]
        public override double MinSequenceMass { get; set; }

        [Option("maxMass", /*Min = 1,*/ HelpText = "Maximum sequence mass in Da")]
        public override double MaxSequenceMass { get; set; }

        [Option("feature", HelpText = "*.ms1ft, *_isos.csv, or *.msalign (Default: Run ProMex)", HelpShowsDefault = false)]
        public override string FeatureFilePath { get; set; }

        [Option("threads", Min = 0, HelpText = "Maximum number of threads, or 0 to set automatically")]
        public override int MaxNumThreads { get; set; }

        [Option("act", "ActivationMethod", HelpText = "Activation Method")]
        public override ActivationMethod ActivationMethod { get; set; }

        [Option("scansFile", HelpText = "Text file with MS2 scans to process", HelpShowsDefault = false)]
        public string ScansFilePath { get; set; }

        [Option("flip", HelpText = "If specified, FLIP scoring code will be used\n(supports UVPD spectra)")]
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
        }

        /// <summary>
        /// Validate options
        /// </summary>
        /// <parameterFilePath>Path to the parameter file specified by /ParamFile or -ParamFile, or an empty string if not provided</parameterFilePath>
        /// <returns>True if options are valid, otherwise false</returns>
        /// <remarks>Will also load dynamic and static modifications from the parameter file, if defined</remarks>
        public bool Validate(string parameterFilePath)
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

            // This will be True if SpecFilePath is a directory that is NOT a supported folder-type dataset.
            var specPathIsDirectory = Directory.Exists(SpecFilePath) && !isDirectoryDataset;

            if (!File.Exists(SpecFilePath) && !specPathIsDirectory && !isDirectoryDataset)
            {
                PrintError("File not found: " + SpecFilePath);
                return false;
            }

            // The extensions in this variable start with a period
            var supportedFileExtensions = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;

            if (!specPathIsDirectory && !(supportedFileExtensions.Select(ext => SpecFilePath.ToLower().EndsWith(ext)).Any()))
            {
                PrintError("Invalid file extension for spectrum file (" + Path.GetExtension(SpecFilePath) + "): " + SpecFilePath);
                return false;
            }

            SpecFilePaths = GetSpecFilePaths(SpecFilePath);

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

            if (!FastaDatabaseConstants.ValidFASTAExtension(DatabaseFilePath))
            {
                PrintError("Invalid extension for the database file path (" + Path.GetExtension(DatabaseFilePath) + ")");
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
                        return false;
                }
            }
            else
            {
                if (!File.Exists(ModsFilePath))
                {
                    PrintError("Modifications file not found: " + ModsFilePath);
                    return false;
                }

                var errorMessage = LoadModsFile(ModsFilePath);
                if (!string.IsNullOrWhiteSpace(errorMessage))
                {
                    PrintError(errorMessage);
                    return false;
                }
            }

            // Scans file validation
            if (!string.IsNullOrWhiteSpace(ScansFilePath) && !File.Exists(ScansFilePath))
            {
                PrintError("Scans File not found: " + ScansFilePath);
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
                !Path.GetExtension(FeatureFilePath).Equals(".csv", StringComparison.OrdinalIgnoreCase) &&
                !Path.GetExtension(FeatureFilePath).Equals(".ms1ft", StringComparison.OrdinalIgnoreCase) &&
                !Path.GetExtension(FeatureFilePath).Equals(".msalign", StringComparison.OrdinalIgnoreCase))
            {
                PrintError("Invalid extension for the Feature file path (" + Path.GetExtension(FeatureFilePath) + ")");
                return false;
            }

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

        private static void PrintError(string errorMessage, Exception ex = null)
        {
            Console.WriteLine();
            ConsoleMsgUtils.ShowWarning(
                "----------------------------------------------------------\n" +
                "Error: " + errorMessage + "\n" +
                "----------------------------------------------------------");

            if (ex != null)
            {
                ConsoleMsgUtils.ShowWarning(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
            }

            Console.WriteLine();
        }

        public void Display(string parameterFilePath)
        {
            Console.WriteLine();

            Console.WriteLine("MaxThreads: " + MaxNumThreads);

            foreach (var specFilePath in SpecFilePaths)
            {
                Console.WriteLine("SpectrumFilePath:  " + specFilePath);
            }

            Console.WriteLine("DatabaseFilePath:  " + DatabaseFilePath);
            Console.WriteLine("DatabaseFilePath:  " + DatabaseFilePath);
            Console.WriteLine("FeatureFilePath:   {0}", FeatureFilePath ?? "N/A");

            var parameterFilePathToShow = string.IsNullOrWhiteSpace(parameterFilePath) ? "N/A" : parameterFilePath;
            Console.WriteLine("ParameterFilePath: {0}", parameterFilePathToShow);

            Console.WriteLine("OutputDir:         " + OutputDir);
            Console.WriteLine();
            Console.WriteLine("InternalCleavageMode:       " + InternalCleavageMode);
            Console.WriteLine("Tag-based search:           " + TagBasedSearch);
            Console.WriteLine("Tda:                        " + (TargetDecoySearchMode == DatabaseSearchMode.Both ? "Target+Decoy" : TargetDecoySearchMode.ToString()));
            Console.WriteLine("PrecursorIonTolerancePpm:   " + PrecursorIonTolerancePpm);
            Console.WriteLine("ProductIonTolerancePpm:     " + ProductIonTolerancePpm);
            Console.WriteLine("MinSequenceLength:          " + MinSequenceLength);
            Console.WriteLine("MaxSequenceLength:          " + MaxSequenceLength);
            Console.WriteLine("MinPrecursorIonCharge:      " + MinPrecursorIonCharge);
            Console.WriteLine("MaxPrecursorIonCharge:      " + MaxPrecursorIonCharge);
            Console.WriteLine("MinProductIonCharge:        " + MinProductIonCharge);
            Console.WriteLine("MaxProductIonCharge:        " + MaxProductIonCharge);
            Console.WriteLine("MinSequenceMass:            " + MinSequenceMass);
            Console.WriteLine("MaxSequenceMass:            " + MaxSequenceMass);
            Console.WriteLine("MatchesPerSpectrumToReport: " + MatchesPerSpectrumToReport);
            Console.WriteLine("IncludeDecoyResults:        " + IncludeDecoyResults);
            Console.WriteLine("MaxDynamicModificationsPerSequence: " + MaxDynamicModificationsPerSequence);
            Console.WriteLine("Modifications:");

            foreach (var searchMod in Modifications)
            {
                Console.WriteLine(searchMod.ToString(true));
            }

            if (!string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                Console.WriteLine("Getting MS1 features from " + FeatureFilePath);
            }

            if (ScanNumbers?.Any() == true)
            {
                Console.WriteLine("Processing specific MS2 scans:");
                Console.WriteLine(string.Join(", ", ScanNumbers));
            }

            if (UseFLIP)
                Console.WriteLine("Using FLIP scoring.");

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

        private static IEnumerable<string> GetSpecFilePaths(string fileOrDirectoryPath)
        {
            try
            {

                if (Directory.Exists(fileOrDirectoryPath) &&
                    !MassSpecDataReaderFactory.IsADirectoryDataset(fileOrDirectoryPath))
                {
                    // Process all .raw files in the given directory
                    return Directory.GetFiles(fileOrDirectoryPath, "*.raw");
                }
            }
            catch (Exception ex)
            {
                PrintError(string.Format("Exception examining the file or directory path ({0}): {1}", fileOrDirectoryPath, ex.Message), ex);
            }

            // Use the file or directory path as-is
            return new[] { fileOrDirectoryPath };
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
                    return "Error while parsing " + modFilePath;

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
                return string.Empty;

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
                        continue;

                    var dataValues = dataLine.Split(delimiters);
                    foreach (var value in dataValues)
                    {
                        if (!int.TryParse(value, out var scanNumber))
                            continue;

                        if (!scanNumbers.Contains(scanNumber))
                            scanNumbers.Add(scanNumber);
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
                            var errMsg = "Invalid value for NumMods in the MSPathFinder parameter file";
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
                PrintError("Exception extracting dynamic and static mod information from the MSPathFinder parameter file", ex);
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
                            var errMsg = "Static mod definition contains ',opt,'; update the param file to have ',fix,' or change to 'DynamicMod='";
                            PrintError(errMsg + "\n" + staticMod);
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
                            var errMsg = "Dynamic mod definition contains ',fix,'; update the param file to have ',opt,' or change to 'StaticMod='";
                            PrintError(errMsg + "\n" + dynamicMod);
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
                PrintError("Exception storing dynamic and static mod information loaded from the MSPathFinder parameter file", ex);
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
            PrintError(errorMessage, ex);
        }

        private static void WarningEventHandler(string warningMessage)
        {
            ConsoleMsgUtils.ShowWarning(warningMessage);
        }

        #endregion
    }
}
