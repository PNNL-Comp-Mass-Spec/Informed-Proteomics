using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;
using PRISM;

namespace PbfGen
{
    public class PbfGenInputParameters
    {
        [Option("i", "s", ArgPosition = 1, Required = true, HelpShowsDefault = false,
            HelpText = "Input file path: .raw or a directory with .raw files (also supports other input formats, see documentation)")]
        public string SourcePath { get; set; }

        [Option("o", HelpShowsDefault = false,
            HelpText = "Output directory. Defaults to directory containing input file.")]
        public string OutputDir { get; set; }

        [Option("start",
            HelpText = "Start scan number (use to limit scan range included in .pbf file).", Min = -1)]
        public int StartScan { get; set; }

        [Option("end",
            HelpText = "End scan number (use to limit scan range included in .pbf file).", Min = -1)]
        public int EndScan { get; set; }

        /// <summary>
        /// Source dataset file or directory paths
        /// </summary>
        /// <remarks>
        /// This list is populated by the call to Validate
        /// </remarks>
        public List<FileSystemInfo> SourceDatasetPaths { get;  }

        /// <summary>
        /// Constructor
        /// </summary>
        public PbfGenInputParameters()
        {
            SourcePath = string.Empty;
            OutputDir = string.Empty;
            StartScan = -1;
            EndScan = -1;

            SourceDatasetPaths = new List<FileSystemInfo>();
        }

        /// <summary>
        /// Validate the processing options
        /// </summary>
        /// <remarks>This method will populate SourceDatasetPaths</remarks>
        /// <returns>True if success, false if an error</returns>
        public bool Validate()
        {
            var defaultOutputDirectoryPath = string.Empty;
            SourceDatasetPaths.Clear();

            if (SourcePath.Contains("*") || SourcePath.Contains("?"))
            {
                // SourcePath has wildcards
                // Validate each matching file or directory

                var cleanPath = SourcePath.Replace('*', '_').Replace('?','_');
                var lastSlash = SourcePath.LastIndexOf(Path.DirectorySeparatorChar);

                string wildcardSpec;
                if (lastSlash >= 0)
                {
                    wildcardSpec = SourcePath.Substring(lastSlash + 1);
                }
                else
                {
                    wildcardSpec = SourcePath;
                }

                var sourcePathInfo = new FileInfo(cleanPath);
                var workingDirectory = sourcePathInfo.Directory ?? new DirectoryInfo(".");

                var matchingFiles = workingDirectory.GetFiles(wildcardSpec).ToList();
                if (matchingFiles.Count > 0)
                {
                    Console.WriteLine("Finding dataset files that match " + SourcePath);
                    foreach (var datasetFile in matchingFiles)
                    {
                        ShowDebug(PathUtils.CompactPathString(datasetFile.FullName, 80), 0);
                        var validDataset = ValidateSourceData(datasetFile.FullName, out _);
                        if (!validDataset)
                            return false;

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
                        Console.WriteLine("Finding dataset directories that match " + SourcePath);
                        foreach (var datasetDirectory in matchingDirectories)
                        {
                            ShowDebug(PathUtils.CompactPathString(datasetDirectory.FullName, 80), 0);
                            var validDataset = ValidateSourceData(datasetDirectory.FullName, out _);
                            if (!validDataset)
                                return false;

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
                var validDataset = ValidateSourceData(SourcePath, out var specPathIsDirectory);
                if (!validDataset)
                    return false;

                Console.WriteLine("{0,-18} {1}", "Input Dataset:", SourcePath);

                if (specPathIsDirectory)
                {
                    var datasetDirectory = new DirectoryInfo(SourcePath);
                    SourceDatasetPaths.Add(datasetDirectory);
                    defaultOutputDirectoryPath = datasetDirectory.FullName;
                }
                else
                {
                    var datasetFile = new FileInfo(SourcePath);
                    SourceDatasetPaths.Add(datasetFile);
                    defaultOutputDirectoryPath = datasetFile.DirectoryName ?? ".";
                }
            }

            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                OutputDir = defaultOutputDirectoryPath;
            }

            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                ShowError("Invalid output file directory: " + OutputDir);
                return false;
            }

            Console.WriteLine("{0,-18} {1}", "Output Directory:", OutputDir);

            if (EndScan < StartScan && EndScan > -1)
            {
                ShowError($"End scan ({EndScan}) must not be less than start scan ({StartScan})!");
                return false;
            }

            if (!Directory.Exists(OutputDir))
            {
                if (File.Exists(OutputDir) && !File.GetAttributes(OutputDir).HasFlag(FileAttributes.Directory))
                {
                    ShowError("OutputDir \"" + OutputDir + "\" is not a directory!");
                    return false;
                }

                ShowDebug("Creating missing output directory");
                Directory.CreateDirectory(OutputDir);
            }

            return true;
        }

        /// <summary>
        /// Validate that the dataset is a valid mass spec dataset
        /// </summary>
        /// <param name="sourceDatasetPath">Path to a dataset file or directory</param>
        /// <param name="specPathIsDirectory">Output: true if the specified path is a directory</param>
        /// <returns>True if valid, false if an error</returns>
        private bool ValidateSourceData(string sourceDatasetPath, out bool specPathIsDirectory)
        {
            // Check for folder-type datasets, and replace specFilePath with the directory name if it is.
            sourceDatasetPath = MassSpecDataReaderFactory.GetDatasetName(sourceDatasetPath);

            var isDirectoryDataset = MassSpecDataReaderFactory.IsADirectoryDataset(sourceDatasetPath);

            // True if sourceDatasetPath is a directory that is NOT a supported folder-type dataset.
            specPathIsDirectory = Directory.Exists(sourceDatasetPath) && !isDirectoryDataset;

            if (!File.Exists(sourceDatasetPath) && !specPathIsDirectory && !isDirectoryDataset)
            {
                ShowError("File not found: " + sourceDatasetPath);
                return false;
            }

            // The extensions in this variable start with a period
            var supportedFileExtensions = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;
            supportedFileExtensions.Remove(".pbf");

            if (!specPathIsDirectory && !supportedFileExtensions.Select(ext => sourceDatasetPath.ToLower().EndsWith(ext)).Any())
            {
                ShowError("Invalid file extension for spectrum file (" + Path.GetExtension(sourceDatasetPath) + "): " + sourceDatasetPath);
                return false;
            }

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
