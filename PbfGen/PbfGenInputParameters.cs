using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;
using PRISM;

namespace PbfGen
{
    public class PbfGenInputParameters
    {
        [Option("s", ArgPosition = 1, Required = true, HelpText = "Raw file path: *.raw or directory. (also supports other input formats, see documentation)", HelpShowsDefault = false)]
        public string SourcePath { get; set; }

        [Option("o", HelpText = "Output directory. Defaults to directory containing input file.", HelpShowsDefault = false)]
        public string OutputDir { get; set; }

        [Option("start", HelpText = "Start scan number (use to limit scan range included in .pbf file).", Min = -1)]
        public int StartScan { get; set; }

        [Option("end", HelpText = "End scan number (use to limit scan range included in .pbf file).", Min = -1)]
        public int EndScan { get; set; }

        public PbfGenInputParameters()
        {
            SourcePath = string.Empty;
            OutputDir = string.Empty;
            StartScan = -1;
            EndScan = -1;
        }

        public bool Validate()
        {
            // Check for folder-type datasets, and replace specFilePath with the directory name if it is.
            SourcePath = MassSpecDataReaderFactory.GetDatasetName(SourcePath);

            var isDirectoryDataset = MassSpecDataReaderFactory.IsADirectoryDataset(SourcePath);
            // True if specFilePath is a directory that is NOT a supported folder-type dataset.
            var specPathIsDirectory = Directory.Exists(SourcePath) && !isDirectoryDataset;

            if (!File.Exists(SourcePath) && !specPathIsDirectory && !isDirectoryDataset)
            {
                PrintError("File not found: " + SourcePath);
                return false;
            }

            var types = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;
            types.Remove(".pbf");

            if (!specPathIsDirectory && !types.Select(ext => SourcePath.ToLower().EndsWith(ext)).Any())
            {
                PrintError("Invalid file extension: (" + Path.GetExtension(SourcePath) + ") " + SourcePath);
                return false;
            }

            if (string.IsNullOrWhiteSpace(OutputDir))
            {
                // Must use "Path.GetFullPath" to return the absolute path when the source file is in the working directory
                // But, it could cause problems with too-long paths.
                OutputDir = specPathIsDirectory ? SourcePath : Path.GetDirectoryName(Path.GetFullPath(SourcePath));
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
                Directory.CreateDirectory(OutputDir);
            }

            return true;
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
