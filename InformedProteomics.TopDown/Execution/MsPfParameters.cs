using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Execution
{
    /// <summary>
    /// A class representing a parsed MSPathFinder parameter file.
    /// </summary>
    public class MsPfParameters
    {
        // Ignore Spelling: Tda, tri

        public const string ParameterFileExtension = ".param";

        /// <summary>
        /// Initializes a new instance of the <see cref="MsPfParameters"/> class, with the parameters specifying required options for a search
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="dbFilePath"></param>
        /// <param name="outputDir"></param>
        /// <param name="aaSet"></param>
        /// <param name="featureFilePath"></param>
        public MsPfParameters(string specFilePath, string dbFilePath, string outputDir, AminoAcidSet aaSet, string featureFilePath = null) : this()
        {
            // ReSharper disable VirtualMemberCallInConstructor
            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            OutputDir = outputDir;
            FeatureFilePath = featureFilePath;
            // ReSharper restore VirtualMemberCallInConstructor
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="MsPfParameters"/> class.
        /// </summary>
        public MsPfParameters()
        {
            Modifications = new List<SearchModification>();
            SetDefaults();
        }

        private void SetDefaults()
        {
            FeatureFilePath = null;
            MinSequenceLength = 21;
            MaxSequenceLength = 500;
            MinPrecursorIonCharge = 2;
            MaxPrecursorIonCharge = 50;
            MinProductIonCharge = 1;
            MaxProductIonCharge = 20;
            MinSequenceMass = 3000.0;
            MaxSequenceMass = 50000.0;
            PrecursorIonTolerance = new Tolerance(10);
            ProductIonTolerance = new Tolerance(10);
            TargetDecoySearchMode = DatabaseSearchMode.Target;
            InternalCleavageMode = InternalCleavageType.SingleInternalCleavage;

            MaxNumNTermCleavages = 1;
            MaxNumCTermCleavages = 0;
            MaxNumThreads = 0;
            ScanNumbers = null;
            MatchesPerSpectrumToKeepInMemory = 3;
            MatchesPerSpectrumToReport = 1;
            IncludeDecoyResults = false;
            TagBasedSearch = true;

            ActivationMethod = ActivationMethod.Unknown;
            OverwriteExistingResults = false;
        }

        /// <summary>
        /// Gets or sets the output directory for the results
        /// </summary>
        public virtual string OutputDir { get; set; }

        /// <summary>
        /// Gets or sets the AminoAcidSet
        /// </summary>
        public AminoAcidSet AminoAcidSet { get; set; }

        /// <summary>
        /// Maximum number of N-terminal cleavages
        /// </summary>
        /// <remarks>default 1</remarks>
        public int MaxNumNTermCleavages { get; set; }

        /// <summary>
        /// Maximum number of C-terminal cleavages
        /// </summary>
        /// <remarks>default 0</remarks>
        public int MaxNumCTermCleavages { get; set; }

        /// <summary>
        /// Gets or sets the path for the Spec file.
        /// </summary>
        public virtual string SpecFilePath { get; set; }

        /// <summary>
        /// Gets or sets the path for the FASTA file.
        /// </summary>
        public virtual string DatabaseFilePath { get; set; }

        /// <summary>
        /// Gets or sets the path for the MS1 feature file.
        /// </summary>
        public virtual string FeatureFilePath { get; set; }

        /// <summary>
        /// Gets or sets the MSPathFinder search mode.
        /// </summary>
        public virtual InternalCleavageType InternalCleavageMode { get; set; }

        /// <summary>
        /// Gets or sets the DB search mode.
        /// </summary>
        public virtual DatabaseSearchMode TargetDecoySearchMode { get; set; }

        /// <summary>
        /// If false, and existing results are found, the existing results will be used
        /// Set to True to force the target and/or decoy search to be repeated even if _IcTarget.tsv or _IcDecoy.tsv exists
        /// </summary>
        public virtual bool OverwriteExistingResults { get; set; }

        /// <summary>
        /// Gets or sets the precursor ion tolerance.
        /// </summary>
        public Tolerance PrecursorIonTolerance { get; set; }

        /// <summary>
        /// Gets or sets the activation method to use during search.
        /// </summary>
        public virtual ActivationMethod ActivationMethod { get; set; }

        /// <summary>
        /// Gets or sets the precursor ion tolerance in ppm.
        /// </summary>
        public virtual double PrecursorIonTolerancePpm
        {
            get => PrecursorIonTolerance.GetValue();
            set => PrecursorIonTolerance = new Tolerance(value);
        }

        /// <summary>
        /// Gets or sets the product ion tolerance.
        /// </summary>
        public Tolerance ProductIonTolerance { get; set; }

        /// <summary>
        /// Gets or sets the product ion tolerance in ppm.
        /// </summary>
        public virtual double ProductIonTolerancePpm
        {
            get => ProductIonTolerance.GetValue();
            set => ProductIonTolerance = new Tolerance(value);
        }

        /// <summary>
        /// Gets or sets the minimum length of a sequence.
        /// </summary>
        public virtual int MinSequenceLength { get; set; }

        /// <summary>
        /// Gets or sets the maximum length of a sequence.
        /// </summary>
        public virtual int MaxSequenceLength { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible precursor ion charge state.
        /// </summary>
        public virtual int MinPrecursorIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible precursor ion charge state.
        /// </summary>
        public virtual int MaxPrecursorIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible product ion charge state.
        /// </summary>
        public virtual int MinProductIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible product ion charge state.
        /// </summary>
        public virtual int MaxProductIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible sequence mass.
        /// </summary>
        public virtual double MinSequenceMass { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible sequence mass.
        /// </summary>
        public virtual double MaxSequenceMass { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible MS1 feature probability threshold.
        /// </summary>
        public double MinFeatureProbability { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible modification combinations per sequence.
        /// </summary>
        public int MaxDynamicModificationsPerSequence { get; set; }

        /// <summary>
        /// Gets or sets the DB search mode, using a tri-state boolean
        /// </summary>
        /// <remarks>default true
        /// true: target and decoy, false: target only, null: decoy only</remarks>
        [Obsolete("Use TargetDecoySearchMode")]
        public bool? RunTargetDecoyAnalysisBool
        {
            get
            {
                if (TargetDecoySearchMode == DatabaseSearchMode.Both)
                {
                    return true;
                }

                if (TargetDecoySearchMode == DatabaseSearchMode.Decoy)
                {
                    return null;
                }
                //(Tda2 == DatabaseSearchMode.Target)
                return false;
            }
            set
            {
                if (value == null)
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Decoy;
                }
                else if (value.Value)
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Both;
                }
                else
                {
                    TargetDecoySearchMode = DatabaseSearchMode.Target;
                }
            }
        }

        /// <summary>
        /// Gets or sets the DB search mode
        /// </summary>
        [Obsolete("Use TargetDecoySearchMode")]
        public DatabaseSearchMode RunTargetDecoyAnalysis
        {
            get => TargetDecoySearchMode;
            set => TargetDecoySearchMode = value;
        }

        public virtual bool TagBasedSearch { get; set; }

        /// <summary>
        /// Specific MS2 scan numbers to process
        /// </summary>
        public IEnumerable<int> ScanNumbers { get; set; }

        /// <summary>
        /// Number of matches to track in memory for each spectrum
        /// The matches are used when computing confidence scores
        /// </summary>
        /// <remarks>Defaults to 3</remarks>
        [Obsolete("Superseded by MatchesPerSpectrumToKeepInMemory")]
        public int NumMatchesPerSpectrum => MatchesPerSpectrumToKeepInMemory;

        /// <summary>
        /// Number of matches to track in memory for each spectrum
        /// The matches in memory are used when computing spectral E-values
        /// </summary>
        /// <remarks>Defaults to 3</remarks>
        public virtual int MatchesPerSpectrumToKeepInMemory { get; set; }

        /// <summary>
        /// Number of matches per spectrum to list in the results file
        /// </summary>
        /// <remarks>Defaults to 1</remarks>
        public virtual int MatchesPerSpectrumToReport { get; set; }

        /// <summary>
        /// Gets or sets the option to include decoy results in the _IcTda.tsv file
        /// </summary>
        public virtual bool IncludeDecoyResults { get; set; }

        /// <summary>
        /// Maximum number of threads
        /// </summary>
        /// <remarks>Defaults to 4</remarks>
        public virtual int MaxNumThreads { get; set; }

        /// <summary>
        /// Maximum degree of parallelism (MaxDOP), i.e.
        /// the maximum number of concurrent tasks when running a parallel operation
        /// </summary>
        /// <remarks>
        /// Equivalent to MaxNumThreads if MaxNumThreads is greater than 0.
        /// Otherwise, -1, meaning to let the system decide the number of threads to use
        /// </remarks>
        public int MaxDegreeOfParallelism => MaxNumThreads > 0 ? MaxNumThreads : -1;

        /// <summary>
        /// 0: all internal sequences,
        /// 1: #NCleavages &lt;= Max OR Cleavages &lt;= Max (Default)
        /// 2: 1: #NCleavages &lt;= Max AND Cleavages &lt;= Max
        /// </summary>
        /// <remarks>default 1</remarks>
        [Obsolete("Use InternalCleavageMode")]
        public virtual int SearchModeInt
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

        /// <summary>
        /// Gets or sets a list containing the post-translational modifications used for database search.
        /// </summary>
        public List<SearchModification> Modifications { get; set; }

        public void Write()
        {
            var outputFilePath = Path.Combine(OutputDir, Path.GetFileNameWithoutExtension(SpecFilePath) + ParameterFileExtension);

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("SpecFile\t" + Path.GetFileName(SpecFilePath));
                writer.WriteLine("DatabaseFile\t" + Path.GetFileName(DatabaseFilePath));
                writer.WriteLine("FeatureFile\t{0}", !string.IsNullOrWhiteSpace(FeatureFilePath) ? Path.GetFileName(FeatureFilePath) : Path.GetFileName(MassSpecDataReaderFactory.ChangeExtension(SpecFilePath, ".ms1ft")));
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
                writer.WriteLine("ActivationMethod\t" + ActivationMethod);
                writer.WriteLine("MaxDynamicModificationsPerSequence\t" + MaxDynamicModificationsPerSequence);
                foreach (var searchMod in Modifications)
                {
                    writer.WriteLine("Modification\t" + searchMod);
                }
            }
        }

        /// <summary>
        /// Opens parameter file from .PARAM file..
        /// </summary>
        /// <param name="filePath">The path of the parameter file.</param>
        /// <returns>
        /// Parsed MSPathFinder parameters. Returns null if file does not exist.
        /// Throws exception if file is not formatted correctly.
        /// </returns>
        public static MsPfParameters Parse(string filePath)
        {
            if (!File.Exists(filePath))
            {
                return null;
            }

            var file = File.ReadAllLines(filePath);

            var param = new MsPfParameters();

            foreach (var line in file)
            {
                var parts = line.Split('\t');
                if (parts.Length < 2)
                {
                    continue;
                }

                switch (parts[0])
                {
                    case "SpecFile":
                        param.SpecFilePath = parts[1];
                        break;
                    case "DatabaseFile":
                        param.DatabaseFilePath = parts[1];
                        break;
                    case "FeatureFile":
                        param.FeatureFilePath = parts[1];
                        break;
                    case "SearchMode":
#pragma warning disable 618
                        param.SearchModeInt = Convert.ToInt32(parts[1]);
#pragma warning restore 618
                        break;
                    case "InternalCleavageMode":
                        param.InternalCleavageMode = (InternalCleavageType)Enum.Parse(typeof(InternalCleavageType), parts[1]);
                        break;
                    case "Tda":
                        var tda = 0;
                        tda += Convert.ToInt32(parts[1].Contains("Target"));
                        tda += Convert.ToInt32(parts[1].Contains("Decoy")) * 2;
                        param.TargetDecoySearchMode = (DatabaseSearchMode)tda;
                        break;
                    case "PrecursorIonTolerancePpm":
                        param.PrecursorIonTolerance = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
                        break;
                    case "ProductIonTolerancePpm":
                        param.ProductIonTolerance = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
                        break;
                    case "MinSequenceLength":
                        param.MinSequenceLength = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxSequenceLength":
                        param.MaxSequenceLength = Convert.ToInt32(parts[1]);
                        break;
                    case "MinPrecursorIonCharge":
                        param.MinPrecursorIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxPrecursorIonCharge":
                        param.MaxPrecursorIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MinProductIonCharge":
                        param.MinProductIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxProductIonCharge":
                        param.MaxProductIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MinSequenceMass":
                        param.MinSequenceMass = Convert.ToDouble(parts[1]);
                        break;
                    case "MaxSequenceMass":
                        param.MaxSequenceMass = Convert.ToDouble(parts[1]);
                        break;
                    case "MinFeatureProbability":
                        param.MinFeatureProbability = Convert.ToDouble(parts[1]);
                        break;
                    case "MaxDynamicModificationsPerSequence":
                        param.MaxDynamicModificationsPerSequence = Convert.ToInt32(parts[1]);
                        break;
                    case "Modification":
                        param.Modifications.AddRange(ModFileParser.ParseModification(parts[1]));
                        break;
                    case "ActivationMethod":
                        param.ActivationMethod = (ActivationMethod)Enum.Parse(typeof(ActivationMethod), parts[1]);
                        break;
                }
            }

            if (param.PrecursorIonTolerance == null && param.ProductIonTolerance != null)
            {
                param.PrecursorIonTolerance = param.ProductIonTolerance;
            }

            if (param.PrecursorIonTolerance != null && param.ProductIonTolerance == null)
            {
                param.ProductIonTolerance = param.PrecursorIonTolerance;
            }

            return param;
        }
    }
}
