using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Execution
{
    public class LcMsFeatureFinderInputParameter
    {
        public const string MINIMUM_MASS = "-minMass";
        public const string MAXIMUM_MASS = "-maxMass";
        public const string MINIMUM_CHARGE = "-minCharge";
        public const string MAXIMUM_CHARGE = "-maxCharge";
        public const string INPUT_FILE_PATH = "-i";
        public const string OUTPUT_FOLDER_PATH = "-o";
        public const string MAXIMUM_THREADS = "-maxThreads";
        public const string INCLUDE_ADDITIONAL_SCORES = "-score";
        public const string SAVE_CSV = "-csv";
        public const string SAVE_PNG_FEATURE_MAP = "-featureMap";
        public const string LIKELIHOOD_SCORE_THRESHOLD = "-scoreTh";
        public const string EXISTING_MS1FT_FILE = "-ms1ft";

        public double MinSearchMass;
        public double MaxSearchMass;
        public int MinSearchCharge;
        public int MaxSearchCharge;
        public string InputPath;
        public bool ScoreReport;
        public bool CsvOutput;

        public bool FeatureMapImage;

        public int MaxThreads;
        public string OutputPath;
        public string ExistingFeaturesFilePath;

        public double LikelihoodScoreThreshold;

        public LcMsFeatureFinderInputParameter()
        {
            MinSearchMass = 600;
            MaxSearchMass = 50000;
            MinSearchCharge = 1;
            MaxSearchCharge = 60;
            ScoreReport = false;
            CsvOutput = true;
            LikelihoodScoreThreshold = 0;
            MaxThreads = 0;
        }

        public LcMsFeatureFinderInputParameter(Dictionary<string, string> paramDic)
        {
            Parse(paramDic);
        }

        public void Parse(Dictionary<string, string> paramDic)
        {
            MinSearchMass = Math.Max(double.Parse(paramDic[MINIMUM_MASS]), 600);
            MaxSearchMass = Math.Min(double.Parse(paramDic[MAXIMUM_MASS]), 100000);

            MinSearchCharge = (int)Math.Max(double.Parse(paramDic[MINIMUM_CHARGE]), 1);
            MaxSearchCharge = (int)Math.Min(double.Parse(paramDic[MAXIMUM_CHARGE]), 60);
            InputPath = paramDic[INPUT_FILE_PATH];
            OutputPath = paramDic[OUTPUT_FOLDER_PATH];

            //OutputPath = (OutputPath == null) ? Path.GetDirectoryName(Path.GetFullPath(InputPath)) : Path.GetDirectoryName(OutputPath);
            MaxThreads = int.Parse(paramDic[MAXIMUM_THREADS]);
            MaxThreads = GetOptimalMaxThreads(MaxThreads);

            ScoreReport = Str2Bool(paramDic[INCLUDE_ADDITIONAL_SCORES]);
            CsvOutput = Str2Bool(paramDic[SAVE_CSV]);

            FeatureMapImage = Str2Bool(paramDic[SAVE_PNG_FEATURE_MAP]);

            LikelihoodScoreThreshold = double.Parse(paramDic[LIKELIHOOD_SCORE_THRESHOLD]);

            ExistingFeaturesFilePath = paramDic[EXISTING_MS1FT_FILE];
        }

        public void Display()
        {
            Console.WriteLine("InputPath    {0}", InputPath);
            Console.WriteLine("OutputPath   {0}", OutputPath);

            if (string.IsNullOrEmpty(ExistingFeaturesFilePath))
            {
                Console.WriteLine("MinMass      {0}", MinSearchMass);
                Console.WriteLine("MaxMass      {0}", MaxSearchMass);
                Console.WriteLine("MinCharge    {0}", MinSearchCharge);
                Console.WriteLine("MaxCharge    {0}", MaxSearchCharge);

                Console.WriteLine("FeatureMap   {0}", FeatureMapImage ? "Y" : "N");

                Console.WriteLine("ScoreReport  {0}", ScoreReport ? "Y" : "N");

                Console.WriteLine("LikelihoodRatioThreshold {0}", LikelihoodScoreThreshold);

                Console.WriteLine("MaxThreads   {0}", MaxThreads);
            }
            else
            {
                Console.WriteLine("MS1 features file {0}", ExistingFeaturesFilePath);
            }

        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }

        private int GetOptimalMaxThreads(int userMaxThreads)
        {
            var threads = userMaxThreads;

            if (threads <= 0 || threads > ParallelizationUtils.NumPhysicalCores)
            {
                threads = ParallelizationUtils.NumPhysicalCores;
            }
            return threads;
        }
    }
}
