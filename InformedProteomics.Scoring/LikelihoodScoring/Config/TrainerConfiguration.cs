using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;

namespace InformedProteomics.Scoring.LikelihoodScoring.Config
{
    public enum AcquisitionMethod
    {
        Dia,
        Dda
    };

    public class TrainerConfiguration
    {
        public string[] DataSets { get; private set; }
        public string TsvPath { get; private set; }
        public string DataPath { get; private set; }
        public string OutputPath { get; private set; }
        public DataFileFormat DataFormat { get; private set; }
        
        public int MaxRanks { get; private set; }
        public double RelativeIntensityThreshold { get; private set; }

        public string OutputFileName { get; private set; }

        public IonType[] IonTypes;
        public double SelectedIonThreshold { get; private set; }
        public ActivationMethod ActivationMethod { get; private set; }
        public AcquisitionMethod AcquisitionMethod { get; private set; }
        public int PrecursorCharge { get; private set; }

        public Tolerance Tolerance { get; private set; }
        public Tolerance MassErrorTolerance { get; private set; }

        public int RetentionCount { get; private set; }
        public double WindowWidth { get; private set; }
        public double PrecursorOffsetWidth { get; private set; }
        public double PrecursorOffsetThreshold { get; private set; }
        
        public const double BinWidth = 1.005;
        public int MassBinSize { get; private set; }
        public const double MassErrorBinWidth = 0.01;
        public const double MassErrorWidth = 0.5;
        public const double MassErrorStartPoint = 0.005;

        public int[] SmoothingRanks { get; private set; }
        public int[] SmoothingWindowSize { get; private set; }

        public TrainerConfiguration(string configurationFile)
        {
            ReadConfigurationFile(configurationFile);
        }

        public void ReadConfigurationFile(string configurationFile)
        {
            var reader = new ConfigFileReader(configurationFile);
            // Read program variables
            var config = reader.GetNodes("vars").First();
            PrecursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            PrecursorOffsetThreshold = Convert.ToDouble(config.Contents["precursoroffsetthreshold"]);
            WindowWidth = Convert.ToInt32(config.Contents["searchwidth"]);
            PrecursorOffsetWidth = Convert.ToInt32(config.Contents["precursoroffsetwidth"]);
            RetentionCount = Convert.ToInt32(config.Contents["retentioncount"]);
            RelativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);
            SelectedIonThreshold = Convert.ToDouble(config.Contents["selectedionthreshold"]);
            MassBinSize = Convert.ToInt32(config.Contents["massbinsize"]);
            var actStr = config.Contents["activationmethod"].ToLower();
            switch (actStr)
            {
                case "hcd":
                    ActivationMethod = ActivationMethod.HCD;
                    Tolerance = _defaultTolerancePpm;
                    break;
                case "cid":
                    ActivationMethod = ActivationMethod.CID;
                    Tolerance = _defaultToleranceTh;
                    break;
                case "etd":
                    ActivationMethod = ActivationMethod.ETD;
                    Tolerance = _defaultTolerancePpm;
                    break;
                default:
                    throw new FormatException("Invalid Activation Method.");
            }

            var acqStr = config.Contents["acquisitionmethod"].ToLower();
            switch (acqStr)
            {
                case "dia":
                    AcquisitionMethod = AcquisitionMethod.Dia;
                    break;
                case "dda":
                    AcquisitionMethod = AcquisitionMethod.Dda;
                    break;
                default:
                    throw new FormatException("Invalid Acquisition Method.");
            }

            MassErrorTolerance = _defaultToleranceTh;

            MaxRanks = Convert.ToInt32(config.Contents["maxranks"]);

            var smoothingRanksStr = config.Contents["smoothingranks"].Split(',');
            SmoothingRanks = new int[smoothingRanksStr.Length];
            var smoothingWindowSizeStr = config.Contents["smoothingwindowsize"].Split(',');
            SmoothingWindowSize = new int[smoothingWindowSizeStr.Length];
            if (SmoothingRanks.Length != SmoothingWindowSize.Length)
                throw new ArgumentException("SmoothingRanks and SmoothingWindowSize unequal lengths.");
            for (int i = 0; i < SmoothingRanks.Length; i++)
            {
                if (smoothingRanksStr[i] == "Max") SmoothingRanks[i] = Int32.MaxValue;
                else SmoothingRanks[i] = Convert.ToInt32(smoothingRanksStr[i]);
                SmoothingWindowSize[i] = Convert.ToInt32(smoothingWindowSizeStr[i]);
            }

            // Read ion data
            var ionInfo = reader.GetNodes("ion").First();
            int totalCharges = Convert.ToInt32(ionInfo.Contents["totalcharges"]);
            var ionTypeStr = ionInfo.Contents["iontype"].Split(',');
            var ions = new BaseIonType[ionTypeStr.Length];
            for (int i = 0; i < ionTypeStr.Length; i++)
            {
                switch (ionTypeStr[i].ToLower())
                {
                    case "a":
                        ions[i] = BaseIonType.A;
                        break;
                    case "b":
                        ions[i] = BaseIonType.B;
                        break;
                    case "c":
                        ions[i] = BaseIonType.C;
                        break;
                    case "x":
                        ions[i] = BaseIonType.X;
                        break;
                    case "y":
                        ions[i] = BaseIonType.Y;
                        break;
                    case "z":
                        ions[i] = BaseIonType.Z;
                        break;
                }
            }
            var ionLossStr = ionInfo.Contents["losses"].Split(',');
            var ionLosses = new NeutralLoss[ionLossStr.Length];
            for (int i = 0; i < ionLossStr.Length; i++)
            {
                switch (ionLossStr[i].ToLower())
                {
                    case "noloss":
                        ionLosses[i] = NeutralLoss.NoLoss;
                        break;
                    case "nh3":
                        ionLosses[i] = NeutralLoss.NH3;
                        break;
                    case "h2o":
                        ionLosses[i] = NeutralLoss.H2O;
                        break;
                }
            }
            _ionTypeFactory = new IonTypeFactory(ions, ionLosses, totalCharges);
            IonTypes = _ionTypeFactory.GetAllKnownIonTypes().ToArray();
            var tempIonList = new List<IonType>();
            if (ionInfo.Contents.ContainsKey("exclusions"))
            {
                var ionExclusions = ionInfo.Contents["exclusions"].Split(',');
                tempIonList.AddRange(IonTypes.Where(ionType => !ionExclusions.Contains(ionType.Name)));
                IonTypes = tempIonList.ToArray();
            }

            // Read input and output file names
            var fileInfo = reader.GetNodes("fileinfo").First();
            DataSets = fileInfo.Contents["name"].Split(',');
            var dataFormat = fileInfo.Contents["format"];
            switch (dataFormat)
            {
                case "mgf":
                    DataFormat = DataFileFormat.Mgf;
                    break;
                case "ictopdown":
                    DataFormat = DataFileFormat.IcTopDown;
                    break;
                case "icbottomup":
                    DataFormat = DataFileFormat.IcBottomUp;
                    break;
                case "dia":
                    DataFormat = DataFileFormat.Dia;
                    break;
                default:
                    throw new FormatException("Invalid Acquisition Method.");
            }

            TsvPath = fileInfo.Contents["tsvpath"];
            DataPath = fileInfo.Contents["datapath"];
            var outPathtemp = fileInfo.Contents["outpath"];
            OutputPath = outPathtemp;

            OutputFileName = OutputPath + fileInfo.Contents["outputfile"];
        }

        private readonly Tolerance _defaultTolerancePpm = new Tolerance(10, ToleranceUnit.Ppm);
        private readonly Tolerance _defaultToleranceTh = new Tolerance(0.5, ToleranceUnit.Th);
        private IonTypeFactory _ionTypeFactory;
    }
}
