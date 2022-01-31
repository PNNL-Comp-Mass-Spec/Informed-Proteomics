using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.MassSpecData;
using PRISM;

namespace PbfGen
{
    public static class Program
    {
        // Ignore Spelling: pbf

        public const string Name = "PbfGen";

        public static string Version
        {
            get
            {
                var programVersion = System.Reflection.Assembly.GetExecutingAssembly().GetName().Version;
                return string.Format("version {0}.{1}.{2} (" + InformedProteomics.Backend.Utils.Misc.GetBuildDateTextFromVersion() + ")", programVersion.Major, programVersion.Minor, programVersion.Build);
            }
        }

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        public static int Main(string[] args)
        {
            try
            {
                // Example text:  PbfGen version 1.1.7706 (February 5, 2021)
                // (the build date is computed automatically)

                Console.WriteLine(Name + " " + Version);
                Console.WriteLine();

                var osVersionInfo = new OSVersionInfo();
                if (osVersionInfo.GetOSVersion().IndexOf("windows", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    var handle = Process.GetCurrentProcess().MainWindowHandle;
                    SetConsoleMode(handle, EnableExtendedFlags);
                }

                string entryAsmName;
                try
                {
                    // ReSharper disable once PossibleNullReferenceException
                    entryAsmName = System.Reflection.Assembly.GetEntryAssembly().GetName().Name;
                }
                catch
                {
                    // This method was likely invoked by NUnit
                    entryAsmName = "Unknown_Assembly";
                }

                var parser = new CommandLineParser<PbfGenInputParameters>(entryAsmName, Version);
                parser.UsageExamples.Add($"Using -start and -end to limit the scan range to include in the .pbf file\n\t{entryAsmName}.exe -s Dataset.raw -start 2000 -end 3000");

                var results = parser.ParseArgs(args);

                if (!results.Success)
                {
                    System.Threading.Thread.Sleep(1500);
                    return -1;
                }

                if (!results.ParsedResults.Validate())
                {
                    parser.PrintHelp();

                    System.Threading.Thread.Sleep(1500);
                    return -1;
                }

                var options = results.ParsedResults;

                var errorCode = ProcessFiles(options);

                if (errorCode != 0)
                    System.Threading.Thread.Sleep(1500);

                return errorCode;
            }
            catch (Exception ex)
            {
                ConsoleMsgUtils.ShowError("Exception while parsing the command line parameters", ex);
                System.Threading.Thread.Sleep(1500);
                return -5;
            }
        }

        private static int ProcessFiles(PbfGenInputParameters options)
        {
#if (!DISABLE_TRYCATCH)
            try
#endif
            {
                if (options.SourceDatasetPaths.Count == 0)
                {
                    ConsoleMsgUtils.ShowWarning("Source dataset path(s) could not be determined; nothing to do");
                    return -6;
                }

                var specFilePaths = new List<string>();
                var scanRangeShown = false;

                if (options.SourceDatasetPaths[0] is DirectoryInfo { Exists: true } inputDirectory)
                {
                    if (MassSpecDataReaderFactory.IsADirectoryDataset(inputDirectory.FullName))
                    {
                        ConsoleMsgUtils.ShowDebug("Processing directory based dataset: " + inputDirectory.FullName);
                    }
                    else
                    {
                        // Look for .raw and .mzML files
                        ConsoleMsgUtils.ShowDebug("Looking for .raw and .mzML files in directory " + inputDirectory.FullName);

                        specFilePaths.AddRange(Directory.GetFiles(inputDirectory.FullName, "*.raw"));
                        specFilePaths.AddRange(Directory.GetFiles(inputDirectory.FullName, "*.mzML"));

                        if (specFilePaths.Count == 0)
                        {
                            ConsoleMsgUtils.ShowWarning(
                                "Directory {0} does not have any .raw or .mzML files, and it is not a directory based dataset; aborting",
                                inputDirectory.FullName);

                            return -7;
                        }
                    }
                }

                if (specFilePaths.Count == 0)
                {
                    specFilePaths.AddRange(options.SourceDatasetPaths.Select(item => item.FullName));
                }

                foreach (var datasetFilePath in specFilePaths)
                {
                    Console.WriteLine();

                    var pbfFileName = MassSpecDataReaderFactory.ChangeExtension(datasetFilePath, PbfLcMsRun.FileExtensionConst);

                    if (string.IsNullOrWhiteSpace(pbfFileName))
                    {
                        ConsoleMsgUtils.ShowWarning("MassSpecDataReaderFactory.ChangeExtension could not determine the pbf filename for {0}", datasetFilePath);
                        continue;
                    }

                    var pbfFilePath = Path.Combine(options.OutputDir, Path.GetFileName(pbfFileName));

                    if (File.Exists(pbfFilePath) && PbfLcMsRun.CheckFileFormatVersion(pbfFilePath, out var isCurrent) && isCurrent)
                    {
                        Console.WriteLine("{0} already exists; will not re-generate", pbfFilePath);
                        continue;
                    }

                    Console.WriteLine("Creating {0}", pbfFilePath);
                    Console.WriteLine("from     {0}", datasetFilePath);

                    if (!scanRangeShown)
                    {
                        ShowScanRange(options);
                        scanRangeShown = true;
                    }

                    var reader = MassSpecDataReaderFactory.GetMassSpecDataReader(datasetFilePath);
                    var progress = new Progress<ProgressData>(p =>
                    {
                        p.UpdateFrequencySeconds = 2;
                        if ((p.Percent % 25).Equals(0) || p.ShouldUpdate())
                        {
                            Console.Write("\r{0}, {1:00.0}% complete                        ", p.Status, p.Percent);
                        }
                    });

                    const double precursorSignalToNoiseRatioThreshold = 0.0;
                    const double productSignalToNoiseRatioThreshold = 0.0;
                    const bool keepDataReaderOpen = false;

                    var run = new PbfLcMsRun(
                        datasetFilePath, reader, pbfFilePath,
                        precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold,
                        progress, keepDataReaderOpen, options.StartScan, options.EndScan);

                    Console.WriteLine();
                }

                Console.WriteLine();
                Console.WriteLine("PbfFormatVersion: {0}", PbfLcMsRun.FileFormatId);
                return 0;
            }
#if (!DISABLE_TRYCATCH)
            catch (Exception ex)
            {
                // NOTE: The DMS Analysis Manager looks for this text; do not change it
                ConsoleMsgUtils.ShowError("Exception while processing", ex);

                var errorCode = -Math.Abs(ex.Message.GetHashCode());
                return errorCode == 0 ? -1 : errorCode;
            }
#endif
        }

        private static void ShowScanRange(PbfGenInputParameters options)
        {
            if (options.StartScan > 0 && options.EndScan > 0)
            {
                Console.WriteLine();
                Console.WriteLine("Only including scans {0} to {1}", options.StartScan, options.EndScan);
            }
            else if (options.StartScan > 0)
            {
                Console.WriteLine();
                Console.WriteLine("Only including scans {0} to the end", options.StartScan);
            }
            else if (options.EndScan > 0)
            {
                Console.WriteLine();
                Console.WriteLine("Only including scans 1 to {0}", options.EndScan);
            }
        }
    }
}
