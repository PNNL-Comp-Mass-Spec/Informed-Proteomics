using System;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.FeatureFinding
{
    public class LcMsFeatureFinderInputParameters
    {
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
            Console.WriteLine("InputPath    {0}", InputPath);
            Console.WriteLine("OutputPath   {0}", OutputPath);

            if (string.IsNullOrEmpty(ExistingFeaturesFilePath))
            {
                Console.WriteLine("MinMass     {0,7:N0} Da", MinSearchMass);
                Console.WriteLine("MaxMass     {0,7:N0} Da", MaxSearchMass);
                Console.WriteLine("MinCharge    {0,2}", MinSearchCharge);
                Console.WriteLine("MaxCharge    {0,2}", MaxSearchCharge);

                Console.WriteLine("FeatureMap   {0}", FeatureMapImage);

                Console.WriteLine("ScoreReport  {0}", ScoreReport);

                Console.WriteLine("LikelihoodRatioThreshold {0}", LikelihoodScoreThreshold);

                Console.WriteLine("MaxThreads   {0}", MaxThreads);
            }
            else
            {
                Console.WriteLine("MS1 features file {0}", ExistingFeaturesFilePath);
            }
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
