using System;
using System.IO;
using InformedProteomics.FeatureFinding;
using PRISM;

namespace ProMex
{
    public class ProMexInputParameters : LcMsFeatureFinderInputParameters
    {
        // Ignore Spelling: pbf, Da, heatmap, csv

        [Option("i", ArgPosition = 1, Required = true, HelpText = "Input file or input folder; supports .pbf, .mzML, and several vendor formats (see documentation)", HelpShowsDefault = false)]
        public override string InputPath { get; set; }

        [Option("o", HelpText = "Output directory. (Default: directory containing input file)", HelpShowsDefault = false)]
        public override string OutputPath { get; set; }

        [Option("minCharge", HelpText = "Minimum charge state", Min = 1, Max = 60)]
        public override int MinSearchCharge { get; set; }

        [Option("maxCharge", HelpText = "Maximum charge state", Min = 1, Max = 60)]
        public override int MaxSearchCharge { get; set; }

        [Option("minMass", HelpText = "Minimum mass in Da", Min = 600, Max = 100000)]
        public override double MinSearchMass { get; set; }

        [Option("maxMass", HelpText = "Maximum mass in Da", Min = 600, Max = 100000)]
        public override double MaxSearchMass { get; set; }

        [Option("featureMap", HelpText = "Output the feature heatmap")]
        public override bool FeatureMapImage { get; set; }

        [Option("score", HelpText = "Output extended scoring information")]
        public override bool ScoreReport { get; set; }

        [Option("maxThreads", HelpText = "Max number of threads to use (Default: 0 (automatically determine the number of threads to use))", HelpShowsDefault = false, Min = 0)]
        public override int MaxThreads { get; set; }

        [Option("csv", HelpText = "Also write feature data to a CSV file")]
        public override bool CsvOutput { get; set; }

        [Option("scoreTh", HelpText = "Likelihood score threshold")]
        public override double LikelihoodScoreThreshold { get; set; }

        [Option("ms1ft", HelpText = "ms1ft format feature file path (use '.' to infer the name from the pbf file)", HelpShowsDefault = false)]
        public override string ExistingFeaturesFilePath { get; set; }

        public ProMexInputParameters() : base()
        {
            //InputPath = "";
            //OutputPath = "";
            //MinSearchCharge = 1;
            //MaxSearchCharge = 60;
            //MinSearchMass = 2000;
            //MaxSearchMass = 50000;
            //FeatureMapImage = true;
            //ScoreReport = false;
            //CsvOutput = false;
            //LikelihoodScoreThreshold = -10;
            //MaxThreads = 0;
            //ExistingFeaturesFilePath = "";
        }

        public bool Validate()
        {
            if (!File.Exists(InputPath) && !Directory.Exists(InputPath))
            {
                PrintError("File not found: " + InputPath);
                return false;
            }

            if (MinSearchMass > MaxSearchMass)
            {
                PrintError("minMass must be less than maxMass!");
                return false;
            }

            if (MinSearchCharge > MaxSearchCharge)
            {
                PrintError("minCharge must be less than maxCharge!");
                return false;
            }

            MaxThreads = GetOptimalMaxThreads(MaxThreads);

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
    }
}
