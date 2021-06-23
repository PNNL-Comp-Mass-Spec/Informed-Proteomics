using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.FeatureFinding;
using PRISM;

namespace ProMex
{
    public class ProMexInputParameters : LcMsFeatureFinderInputParameters
    {
        // Ignore Spelling: csv, Da, heatmap, pbf, wildcards

        [Option("InputFile", "i", "s", ArgPosition = 1, Required = true,
            HelpText = "Input file or input directory; supports .pbf, .mzML, and several vendor formats (see documentation)",
            HelpShowsDefault = false)]
        public override string InputPath { get; set; }

        [Option("OutputDirectory", "o",
            HelpText = "Output directory. (Default: directory containing input file)",
            HelpShowsDefault = false)]
        public override string OutputPath { get; set; }

        [Option("MinCharge",
            HelpText = "Minimum charge state", Min = 1, Max = 60)]
        public override int MinSearchCharge { get; set; }

        [Option("MaxCharge",
            HelpText = "Maximum charge state", Min = 1, Max = 60)]
        public override int MaxSearchCharge { get; set; }

        [Option("MinMass",
            HelpText = "Minimum mass in Da", Min = 600, Max = 100000)]
        public override double MinSearchMass { get; set; }

        [Option("MaxMass",
            HelpText = "Maximum mass in Da", Min = 600, Max = 100000)]
        public override double MaxSearchMass { get; set; }

        [Option("FeatureMap",
            HelpText = "Output the feature heatmap. To disable, use -FeatureMap:false or include 'FeatureMap=False' in a parameter file")]
        public override bool FeatureMapImage { get; set; }

        [Option("Score",
            HelpText = "Output extended scoring information")]
        public override bool ScoreReport { get; set; }

        [Option("MaxThreads", Min = 0,
            HelpText = "Max number of threads to use (Default: 0, meaning automatically determine the number of threads to use)",
            HelpShowsDefault = false)]
        public override int MaxThreads { get; set; }

        [Option("csv",
            HelpText = "Also write feature data to a CSV file")]
        public override bool CsvOutput { get; set; }

        [Option("BinningResolutionPPM", "BinResPPM", Min = 1, Max= 128,
            HelpText = "Binning resolution, in ppm. Valid values are 1, 2, 4, 8, 16, 32, 64, or 128")]
        public int BinningResolutionPPM { get; set; } = 16;

        [Option("ScoreThreshold", "ScoreTh",
            HelpText = "Likelihood score threshold")]
        public override double LikelihoodScoreThreshold { get; set; }

        [Option("ms1ft",
            HelpText = "ms1ft format feature file to create graphics files for (the .pbf file must be included in the same directory); " +
                       "use '.' to infer the name from the pbf file",
            HelpShowsDefault = false)]
        public override string ExistingFeaturesFilePath { get; set; }

        /// <summary>
        /// Source dataset file or directory paths
        /// </summary>
        /// <remarks>
        /// This list is populated by the call to Validate
        /// </remarks>
        public List<FileSystemInfo> SourceDatasetPaths { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ProMexInputParameters()
        {
            //InputPath = string.Empty;
            //OutputPath = string.Empty;
            //MinSearchCharge = 1;
            //MaxSearchCharge = 60;
            //MinSearchMass = 2000;
            //MaxSearchMass = 50000;
            //FeatureMapImage = true;
            //ScoreReport = false;
            //CsvOutput = false;
            //LikelihoodScoreThreshold = -10;
            //MaxThreads = 0;
            //ExistingFeaturesFilePath = string.Empty;

            SourceDatasetPaths = new List<FileSystemInfo>();
        }

        public bool Validate()
        {
            var defaultOutputDirectoryPath = string.Empty;
            SourceDatasetPaths.Clear();

            BitCountForBinning = GetBitCountForPPMResolution(BinningResolutionPPM);

            if (InputPath.Contains("*") || InputPath.Contains("?"))
            {
                // InputPath has wildcards
                // Validate each matching file or directory

                var cleanPath = InputPath.Replace('*', '_').Replace('?', '_');
                var lastSlash = InputPath.LastIndexOf(Path.DirectorySeparatorChar);

                string wildcardSpec;
                if (lastSlash >= 0)
                {
                    wildcardSpec = InputPath.Substring(lastSlash + 1);
                }
                else
                {
                    wildcardSpec = InputPath;
                }

                var sourcePathInfo = new FileInfo(cleanPath);
                var workingDirectory = sourcePathInfo.Directory ?? new DirectoryInfo(".");

                var matchingFiles = workingDirectory.GetFiles(wildcardSpec).ToList();
                if (matchingFiles.Count > 0)
                {
                    Console.WriteLine("Finding dataset files that match " + InputPath);
                    foreach (var datasetFile in matchingFiles)
                    {
                        ShowDebug(PathUtils.CompactPathString(datasetFile.FullName, 80), 0);
                        SourceDatasetPaths.Add(datasetFile);
                        defaultOutputDirectoryPath = workingDirectory.FullName;
                    }

                    Console.WriteLine();
                }
                else
                {
                    var matchingDirectories = workingDirectory.GetDirectories(wildcardSpec).ToList();
                    if (matchingDirectories.Count > 0)
                    {
                        Console.WriteLine("Finding dataset directories that match " + InputPath);
                        foreach (var datasetDirectory in matchingDirectories)
                        {
                            ShowDebug(PathUtils.CompactPathString(datasetDirectory.FullName, 80), 0);
                            SourceDatasetPaths.Add(datasetDirectory);
                            defaultOutputDirectoryPath = workingDirectory.FullName;
                        }
                        Console.WriteLine();
                    }
                    else
                    {
                        ShowWarning(string.Format(
                            "No files or directories matched '{0}' in directory {1}", wildcardSpec, workingDirectory.FullName));

                        return false;
                    }
                }
            }
            else
            {
                if (File.Exists(InputPath))
                {
                    var datasetFile = new FileInfo(InputPath);
                    SourceDatasetPaths.Add(datasetFile);
                    defaultOutputDirectoryPath = datasetFile.DirectoryName ?? ".";
                }
                else if (Directory.Exists(InputPath))
                {
                    var datasetDirectory = new DirectoryInfo(InputPath);
                    SourceDatasetPaths.Add(datasetDirectory);
                    defaultOutputDirectoryPath = datasetDirectory.FullName;
                }
                else
                {
                    ShowError("File not found: " + InputPath);
                    return false;
                }
            }

            if (string.IsNullOrWhiteSpace(OutputPath))
            {
                OutputPath = defaultOutputDirectoryPath;
            }

            if (MinSearchMass > MaxSearchMass)
            {
                ShowError("minMass must be less than maxMass!");
                return false;
            }

            if (MinSearchCharge > MaxSearchCharge)
            {
                ShowError("minCharge must be less than maxCharge!");
                return false;
            }

            MaxThreads = GetOptimalMaxThreads(MaxThreads);

            return true;
        }

        private static void ShowDebug(string message, int emptyLinesBeforeMessage = 1)
        {
            ConsoleMsgUtils.ShowDebugCustom(message, emptyLinesBeforeMessage: emptyLinesBeforeMessage);
        }

        private static void ShowError(string errorMessage)
        {
            ShowErrorOrWarning(errorMessage);
        }

        private static void ShowWarning(string message)
        {
            ShowErrorOrWarning(message, string.Empty);
        }

        private static void ShowErrorOrWarning(string message, string messagePrefix = "Error: ")
        {
            Console.WriteLine();
            ConsoleMsgUtils.ShowWarning("{0}\n{1}{2}\n{0}",
                "----------------------------------------------------------",
                messagePrefix, message);
        }
    }
}
