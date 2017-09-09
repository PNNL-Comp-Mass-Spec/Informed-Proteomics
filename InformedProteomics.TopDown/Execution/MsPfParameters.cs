using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.TopDown.Execution
{
    /// <summary>
    /// A class representing a parsed MSPathFinder parameter file.
    /// </summary>
    public class MsPfParameters
    {
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
            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            OutputDir = outputDir;
            FeatureFilePath = featureFilePath;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="MsPfParameters"/> class.
        /// </summary>
        public MsPfParameters()
        {
            this.Modifications = new List<SearchModification>();
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
            TargetDecoySearchMode = DatabaseSearchMode.Both;
            InternalCleavageMode = InternalCleavageType.SingleInternalCleavage;

            MaxNumNTermCleavages = 1;
            MaxNumCTermCleavages = 0;
            MaxNumThreads = 0;
            ScanNumbers = null;
            NumMatchesPerSpectrum = 3;
            TagBasedSearch = true;

            ActivationMethod = ActivationMethod.Unknown;
        }

        /// <summary>
        /// Gets or sets the output directory for the results
        /// </summary>
        public string OutputDir { get; set; }

        /// <summary>
        /// Gets or sets the AminoAcidSet
        /// </summary>
        public AminoAcidSet AminoAcidSet { get; set; }

        /// <remarks>default 1</remarks>
        public int MaxNumNTermCleavages { get; set; }

        /// <remarks>default 0</remarks>
        public int MaxNumCTermCleavages { get; set; }

        /// <summary>
        /// Gets or sets the path for the Spec file.
        /// </summary>
        public string SpecFilePath { get; set; }

        /// <summary>
        /// Gets or sets the path for the FASTA file.
        /// </summary>
        public string DatabaseFilePath { get; set; }

        /// <summary>
        /// Gets or sets the path for the MS1 feature file.
        /// </summary>
        public string FeatureFilePath { get; set; }

        /// <summary>
        /// Gets or sets the MSPathFinder search mode.
        /// </summary>
        public InternalCleavageType InternalCleavageMode { get; set; }

        /// <summary>
        /// Gets or sets the DB search mode.
        /// </summary>
        public DatabaseSearchMode TargetDecoySearchMode { get; set; }

        /// <summary>
        /// Gets or sets the precursor ion tolerance.
        /// </summary>
        public Tolerance PrecursorIonTolerance { get; set; }

        /// <summary>
        /// Gets or sets the activation method to use during search.
        /// </summary>
        public ActivationMethod ActivationMethod { get; set; }

        /// <summary>
        /// Gets or sets the precursor ion tolerance in ppm.
        /// </summary>
        public double PrecursorIonTolerancePpm
        {
            get { return PrecursorIonTolerance.GetValue(); }
            set { PrecursorIonTolerance = new Tolerance(value); }
        }

        /// <summary>
        /// Gets or sets the product ion tolerance.
        /// </summary>
        public Tolerance ProductIonTolerance { get; set; }

        /// <summary>
        /// Gets or sets the product ion tolerance in ppm.
        /// </summary>
        public double ProductIonTolerancePpm
        {
            get { return ProductIonTolerance.GetValue(); }
            set { ProductIonTolerance = new Tolerance(value); }
        }

        /// <summary>
        /// Gets or sets the minimum length of a sequence.
        /// </summary>
        public int MinSequenceLength { get; set; }

        /// <summary>
        /// Gets or sets the maximum length of a sequence.
        /// </summary>
        public int MaxSequenceLength { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible precursor ion charge state.
        /// </summary>
        public int MinPrecursorIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible precursor ion charge state.
        /// </summary>
        public int MaxPrecursorIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible product ion charge state.
        /// </summary>
        public int MinProductIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible product ion charge state.
        /// </summary>
        public int MaxProductIonCharge { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible sequence mass.
        /// </summary>
        public double MinSequenceMass { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible sequence mass.
        /// </summary>
        public double MaxSequenceMass { get; set; }

        /// <summary>
        /// Gets or sets the minimum possible MS1 feature probability threshold.
        /// </summary>
        public double MinFeatureProbablility { get; set; }

        /// <summary>
        /// Gets or sets the maximum possible modification combinations per sequence.
        /// </summary>
        public int MaxDynamicModificationsPerSequence { get; set; }

        /// <summary>
        /// Gets or sets the DB search mode, using a 'tribool'
        /// </summary>
        /// <remarks>default true
        /// true: target and decoy, false: target only, null: decoy only</remarks>
        [Obsolete("Use TargetDecoySearchMode")]
        public bool? RunTargetDecoyAnalysisBool
        {
            get
            {
                if (TargetDecoySearchMode == DatabaseSearchMode.Both)
                    return true;
                if (TargetDecoySearchMode == DatabaseSearchMode.Decoy)
                    return null;
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
            get { return TargetDecoySearchMode; }
            set { TargetDecoySearchMode = value; }
        }

        public bool TagBasedSearch { get; set; }

        /// <summary>
        /// Specific MS2 scan numbers to process
        /// </summary>
        public IEnumerable<int> ScanNumbers { get; set; }

        /// <remarks>default 3</remarks>
        public int NumMatchesPerSpectrum { get; set; }

        /// <remarks>default 4</remarks>
        public int MaxNumThreads { get; set; }

        /// <summary>
        /// 0: all internal sequences,
        /// 1: #NCleavages &lt;= Max OR Cleavages &lt;= Max (Default)
        /// 2: 1: #NCleavages &lt;= Max AND Cleavages &lt;= Max
        /// </summary>
        /// <remarks>default 1</remarks>
        [Obsolete("Use InternalCleavageMode")]
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

        /// <summary>
        /// Gets or sets a list containing the post-translational modifications used for database search.
        /// </summary>
        public List<SearchModification> Modifications { get; set; }

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
                        param.InternalCleavageMode = (InternalCleavageType)Convert.ToInt32(parts[1]);
                        break;
                    case "Tda":
                        int tda = 0;
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
                        param.MinFeatureProbablility = Convert.ToDouble(parts[1]);
                        break;
                    case "MaxDynamicModificationsPerSequence":
                        param.MaxDynamicModificationsPerSequence = Convert.ToInt32(parts[1]);
                        break;
                    case "Modification":
                        param.Modifications.AddRange(ModFileParser.ParseModification(parts[1]));
                        break;
                    case "ActivationMethod":
                        param.ActivationMethod = (ActivationMethod)Convert.ToInt32(parts[1]);
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
