using System;
using System.Collections.Generic;
using System.IO;

namespace InformedProteomics.TopDown.Execution
{
    public class LcMsFeatureFinderInputParameter
    {
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
            MaxThreads = -1;
        }

        public LcMsFeatureFinderInputParameter(Dictionary<string, string> paramDic)
        {
            Parse(paramDic);
        }

        public void Parse(Dictionary<string, string> paramDic)
        {
            MinSearchMass = Math.Max(double.Parse(paramDic["-minMass"]), 600);
            MaxSearchMass = Math.Min(double.Parse(paramDic["-maxMass"]), 100000);

            MinSearchCharge = (int)Math.Max(double.Parse(paramDic["-minCharge"]), 1);
            MaxSearchCharge = (int)Math.Min(double.Parse(paramDic["-maxCharge"]), 60);
            InputPath = paramDic["-i"];
            OutputPath = paramDic["-o"];
            //OutputPath = (OutputPath == null) ? Path.GetDirectoryName(Path.GetFullPath(InputPath)) : Path.GetDirectoryName(OutputPath);
            MaxThreads = Int32.Parse(paramDic["-maxThreads"]);

            ScoreReport = Str2Bool(paramDic["-score"]);
            CsvOutput = Str2Bool(paramDic["-csv"]);

            FeatureMapImage = Str2Bool(paramDic["-featureMap"]);

            LikelihoodScoreThreshold = double.Parse(paramDic["-scoreTh"]);
        }

        public void Display()
        {
            Console.WriteLine("InputPath\t{0}", InputPath);
            Console.WriteLine("OutputPath\t{0}", OutputPath);

            Console.WriteLine("MinMass\t{0}", MinSearchMass);
            Console.WriteLine("MaxMass\t{0}", MaxSearchMass);
            Console.WriteLine("MinCharge\t{0}", MinSearchCharge);
            Console.WriteLine("MaxCharge\t{0}", MaxSearchCharge);

            Console.WriteLine("FeatureMap\t{0}", FeatureMapImage ? "Y" : "N");

            Console.WriteLine("ScoreReport\t{0}", ScoreReport ? "Y" : "N");

            Console.WriteLine("LikelihoodRatioThreshold\t{0}", LikelihoodScoreThreshold);

            Console.WriteLine("MaxThreads\t{0}", MaxThreads);
        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }
    }
}
