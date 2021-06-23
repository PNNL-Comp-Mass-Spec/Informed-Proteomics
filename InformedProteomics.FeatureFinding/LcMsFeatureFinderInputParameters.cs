using System;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.FeatureFinding
{
    public class LcMsFeatureFinderInputParameters
    {
        // Ignore Spelling: Da

        public const int DEFAULT_BIT_COUNT_FOR_BINNING = 27;

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

        /// <summary>
        /// Bit count for binning (value between 24 and 31); defaults to 27
        /// </summary>
        /// <remarks>
        /// 31: 1 ppm resolution
        /// 30: 2 ppm
        /// 29: 4 ppm
        /// 28: 8 ppm
        /// 27: 16 ppm
        /// 26: 32 ppm
        /// 25: 64 ppm
        /// 24: 128 ppm
        /// </remarks>
        public int BitCountForBinning { get; set; }

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
            BitCountForBinning = DEFAULT_BIT_COUNT_FOR_BINNING;
        }

        public void Display()
        {
            Console.WriteLine("{0,-14} {1}", "Input Path:", InputPath);
            Console.WriteLine("{0,-14} {1}", "Output Path:", OutputPath);

            if (string.IsNullOrEmpty(ExistingFeaturesFilePath))
            {
                Console.WriteLine("{0,-14} {1:N0} Da", "Min Mass:", MinSearchMass);
                Console.WriteLine("{0,-14} {1:N0} Da", "Max  Mass:", MaxSearchMass);
                Console.WriteLine("{0,-14} {1}", "Min Charge:", MinSearchCharge);
                Console.WriteLine("{0,-14} {1}", "Max Charge:", MaxSearchCharge);

                Console.WriteLine("{0,-14} {1}", "Feature Map:", FeatureMapImage);

                Console.WriteLine("{0,-14} {1}", "Score Report:", ScoreReport);

                Console.WriteLine("{0,-14} {1}", "Max Threads:", MaxThreads);

                Console.WriteLine("{0,-22} {1}", "Bit Count For Binning:", BitCountForBinning);

                Console.WriteLine("{0,-22} {1} ppm", "Binning Resolution:", GetPPMResolutionForBitCount(BitCountForBinning));

                Console.WriteLine("{0,-25} {1}", "Likelihood Ratio Threshold:", LikelihoodScoreThreshold);
            }
            else
            {
                Console.WriteLine("{0} {1}", "MS1 features file:", ExistingFeaturesFilePath);
            }

            Console.WriteLine();
        }

        /// <summary>
        /// Get the approximate ppm resolution for the given binning bit count
        /// </summary>
        /// <param name="bitCount"></param>
        /// <returns>ppm resolution</returns>
        protected int GetPPMResolutionForBitCount(int bitCount)
        {
            if (bitCount < 24)
                return 128;

            if (bitCount > 31)
                return 1;

            return (int)Math.Pow(2, 31 - bitCount);
        }

        /// <summary>
        /// Convert ppmResolution to nearest power of 2, then convert to the appropriate bit count
        /// </summary>
        /// <param name="ppmResolution"></param>
        /// <returns>Bit count</returns>
        /// <remarks>
        /// 1 ppm converts to 31 bins
        /// 2 ppm converts to 30 bins
        /// 4 ppm converts to 29 bins
        /// etc.
        /// </remarks>
        public static int GetBitCountForPPMResolution(int ppmResolution)
        {
            if (ppmResolution < 1)
                return DEFAULT_BIT_COUNT_FOR_BINNING;

            if (ppmResolution > 128)
                ppmResolution = 128;

            var closestPPM = 0;
            var closestDelta = int.MaxValue;

            for (var comparisonPPM = 1; comparisonPPM <= 128; comparisonPPM *= 2)
            {
                var delta = Math.Abs(ppmResolution - comparisonPPM);
                if (delta >= closestDelta)
                    continue;

                closestDelta = delta;
                closestPPM = comparisonPPM;
            }

            if (closestPPM == 0)
            {
                closestPPM = 16;
            }

            return 31 - (int)Math.Floor(Math.Log(closestPPM, 2));
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
