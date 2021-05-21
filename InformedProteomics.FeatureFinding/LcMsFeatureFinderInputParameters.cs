using System;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.FeatureFinding
{
    public class LcMsFeatureFinderInputParameters
    {
        // Ignore Spelling: Da

        public virtual double MinSearchMass { get; set; }

        public virtual double MaxSearchMass { get; set; }

        public virtual int MinSearchCharge { get; set; }

        public virtual int MaxSearchCharge { get; set; }

        public virtual string InputPath { get; set; }

        public virtual bool ScoreReport { get; set; }

        public virtual bool CsvOutput { get; set; }

        public virtual bool FeatureMapImage { get; set; }

        public virtual int MaxThreads { get; set; }

        public virtual string OutputPath { get; set; }

        public virtual string ExistingFeaturesFilePath { get; set; }

        public virtual double LikelihoodScoreThreshold { get; set; }

        public LcMsFeatureFinderInputParameters()
        {
            SetDefaults();
        }

        protected void SetDefaults()
        {
            MinSearchMass = 2000;
            MaxSearchMass = 50000;
            MinSearchCharge = 1;
            MaxSearchCharge = 60;
            ScoreReport = false;
            CsvOutput = false;
            LikelihoodScoreThreshold = -10;
            MaxThreads = 0;
            FeatureMapImage = true;
        }

        public void Display()
        {
            Console.WriteLine("{0,-14} {1}", "InputPath:", InputPath);
            Console.WriteLine("{0,-14} {1}", "OutputPath:", OutputPath);

            if (string.IsNullOrEmpty(ExistingFeaturesFilePath))
            {
                Console.WriteLine("{0,-14} {1:N0} Da", "MinMass:", MinSearchMass);
                Console.WriteLine("{0,-14} {1:N0} Da", "MaxMass:", MaxSearchMass);
                Console.WriteLine("{0,-14} {1}", "MinCharge:", MinSearchCharge);
                Console.WriteLine("{0,-14} {1}", "MaxCharge:", MaxSearchCharge);

                Console.WriteLine("{0,-14} {1}", "FeatureMap:", FeatureMapImage);

                Console.WriteLine("{0,-14} {1}", "ScoreReport:", ScoreReport);

                Console.WriteLine("{0,-14} {1}", "MaxThreads:", MaxThreads);

                Console.WriteLine("{0,-14} {1}", "LikelihoodRatioThreshold:", LikelihoodScoreThreshold);
            }
            else
            {
                Console.WriteLine("{0} {1}", "MS1 features file:", ExistingFeaturesFilePath);
            }

            Console.WriteLine();
        }

        protected int GetOptimalMaxThreads(int userMaxThreads)
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
