namespace InformedProteomics.TopDown.Execution
{
    using System;
    using System.Collections.Generic;
    using System.IO;
    using InformedProteomics.Backend.Data.Sequence;
    using InformedProteomics.Backend.Data.Spectrometry;
    using InformedProteomics.Backend.Database;

    /// <summary>
    /// A class representing a parsed MSPathFinder parameter file.
    /// </summary>
    public class MsPfParameters
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="MsPfParameters"/> class.
        /// </summary>
        public MsPfParameters()
        {
            this.Modifications = new List<SearchModification>();
        }

        /// <summary>
        /// Gets or sets the path for the PU Spec file.
        /// </summary>
        public string PuSpecFile { get; set; }

        /// <summary>
        /// Gets or sets the path for the FASTA file.
        /// </summary>
        public string DatabaseFile { get; set; }

        /// <summary>
        /// Gets or sets the path for the MS1 feature file.
        /// </summary>
        public string FeatureFile { get; set; }

        /// <summary>
        /// Gets or sets the MSPathFinder search mode.
        /// </summary>
        public InternalCleavageType SearchMode { get; set; }

        /// <summary>
        /// Gets or sets the TDA.
        /// </summary>
        public DatabaseSearchMode Tda { get; set; }

        /// <summary>
        /// Gets or sets the precursor ion tolerance in ppm.
        /// </summary>
        public Tolerance PrecursorTolerancePpm { get; set; }

        /// <summary>
        /// Gets or sets the product ion tolerance in ppm.
        /// </summary>
        public Tolerance ProductIonTolerancePpm { get; set; }

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
                        param.PuSpecFile = parts[1];
                        break;
                    case "DatabaseFile":
                        param.DatabaseFile = parts[1];
                        break;
                    case "FeatureFile":
                        param.FeatureFile = parts[1];
                        break;
                    case "SearchMode":
                        param.SearchMode = (InternalCleavageType)Convert.ToInt32(parts[1]);
                        break;
                    case "Tda":
                        int tda = 0;
                        tda += Convert.ToInt32(parts[1].Contains("Target"));
                        tda += Convert.ToInt32(parts[1].Contains("Decoy")) * 2;
                        param.Tda = (DatabaseSearchMode)tda;
                        break;
                    case "PrecursorIonTolerancePpm":
                        param.PrecursorTolerancePpm = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
                        break;
                    case "ProductIonTolerancePpm":
                        param.ProductIonTolerancePpm = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
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
                }
            }

            if (param.PrecursorTolerancePpm == null && param.ProductIonTolerancePpm != null)
            {
                param.PrecursorTolerancePpm = param.ProductIonTolerancePpm;
            }

            if (param.PrecursorTolerancePpm != null && param.ProductIonTolerancePpm == null)
            {
                param.ProductIonTolerancePpm = param.PrecursorTolerancePpm;
            }

            return param;
        }
    }
}
