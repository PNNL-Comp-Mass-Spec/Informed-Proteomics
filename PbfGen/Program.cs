using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace PbfGen
{
    public class Program
    {
        public static string Version
        {
            get
            {
                var programVersion = System.Reflection.Assembly.GetExecutingAssembly().GetName().Version;
                return string.Format("version {0}.{1}.{2} (" + Misc.GetBuildDateTextFromVersion() + ")", programVersion.Major, programVersion.Minor, programVersion.Build);
            }
        }

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        public static int Main(string[] args)
        {
            string specFilePath;
            string outputDir;
            var startScan = 0;
            var endScan = 0;

            try
            {
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, EnableExtendedFlags);

                if (args.Length == 0)
                {
                    PrintUsageInfo();
                    return -1;
                }

                var paramDic = new Dictionary<string, string>(StringComparer.InvariantCultureIgnoreCase)
                {
                    {"s", null},
                    {"o", null},
                    {"start", null},
                    {"end", null}
                };

                var arguments = new List<string>();
                if (args.Length % 2 != 0 && !args[0].StartsWith("-"))
                {
                    // Assume the first argument is the spec file path
                    arguments.Add("-s");
                }

                arguments.AddRange(args);

                var switchChars = new char[] { '-', '/' };

                for (var i = 0; i < arguments.Count / 2; i++)
                {
                    var key = arguments[2 * i].TrimStart(switchChars);
                    var value = arguments[2 * i + 1];
                    if (!paramDic.ContainsKey(key))
                    {
                        PrintUsageInfo("Invalid parameter: " + key);
                        return -1;
                    }
                    paramDic[key] = value;
                }

                // Parse command line parameters
                specFilePath = paramDic["s"];

                if (specFilePath == null)
                {
                    PrintUsageInfo("Missing required parameter -s!");
                    return -1;
                }

                // Check for folder-type datasets, and replace specFilePath with the directory name if it is.
                specFilePath = MassSpecDataReaderFactory.GetDatasetName(specFilePath);

                var isDirectoryDataset = MassSpecDataReaderFactory.IsADirectoryDataset(specFilePath);
                // True if specFilePath is a directory that is NOT a supported folder-type dataset.
                var specPathIsDirectory = Directory.Exists(specFilePath) && !isDirectoryDataset;

                if (!File.Exists(specFilePath) && !specPathIsDirectory && !isDirectoryDataset)
                {
                    PrintUsageInfo("File not found: " + specFilePath);
                    return -1;
                }

                var types = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;
                types.Remove(".pbf");

                if (!specPathIsDirectory && !(types.Select(ext => specFilePath.ToLower().EndsWith(ext)).Any()))
                {
                    PrintUsageInfo("Invalid file extension: (" + Path.GetExtension(specFilePath) + ") " + specFilePath);
                    return -1;
                }

                // Must use "Path.GetFullPath" to return the absolute path when the source file is in the working directory
                // But, it could cause problems with too-long paths.
                outputDir = paramDic["o"] ?? (specPathIsDirectory ? specFilePath : Path.GetDirectoryName(Path.GetFullPath(specFilePath)));
                if (outputDir == null)
                {
                    PrintUsageInfo("Invalid output file directory: " + specFilePath);
                    return -1;
                }

                if (paramDic.TryGetValue("start", out var startScanText))
                {
                    if (!int.TryParse(startScanText, out startScan))
                    {
                        PrintUsageInfo("value for -start must be an integer, not " + startScanText);
                        return -1;
                    }
                }

                if (paramDic.TryGetValue("end", out var endScanText))
                {
                    if (!int.TryParse(endScanText, out endScan))
                    {
                        PrintUsageInfo("value for -end must be an integer, not " + endScanText);
                        return -1;
                    }
                }

                if (!Directory.Exists(outputDir))
                {
                    if (File.Exists(outputDir) && !File.GetAttributes(outputDir).HasFlag(FileAttributes.Directory))
                    {
                        PrintUsageInfo("OutputDir " + outputDir + " is not a directory!");
                        return -1;
                    }
                    Directory.CreateDirectory(outputDir);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception while parsing the command line parameters: " + ex.Message);
                return -5;
            }

#if (!DEBUG)
            try
            {
#endif
            var specFilePaths = new[] { specFilePath };
            if (Directory.Exists(specFilePath) && !MassSpecDataReaderFactory.IsADirectoryDataset(specFilePath))
            {
                specFilePaths = Directory.GetFiles(specFilePath, "*.raw"); // TODO: Support folders with other formats in them too...
            }

            foreach (var rawFilePath in specFilePaths)
            {
                var pbfFileName = MassSpecDataReaderFactory.ChangeExtension(rawFilePath, PbfLcMsRun.FileExtensionConst);
                var pbfFilePath = Path.Combine(outputDir, Path.GetFileName(pbfFileName));

                bool isCurrent;
                if (File.Exists(pbfFilePath) && PbfLcMsRun.CheckFileFormatVersion(pbfFilePath, out isCurrent) && isCurrent)
                {
                    Console.WriteLine("{0} already exists.", pbfFilePath);
                    continue;
                }

                Console.WriteLine("Creating {0} from {1}", pbfFilePath, rawFilePath);

                if (startScan > 0 && endScan > 0)
                    Console.WriteLine("Only including scans {0} to {1}", startScan, endScan);
                else if (startScan > 0)
                    Console.WriteLine("Only including scans {0} to the end", startScan);
                else if (endScan > 0)
                    Console.WriteLine("Only including scans 1 to {0}", endScan);

                var reader = MassSpecDataReaderFactory.GetMassSpecDataReader(rawFilePath);
                var progress = new Progress<ProgressData>(p =>
                {
                    p.UpdateFrequencySeconds = 2;
                    if ((p.Percent % 25).Equals(0) || p.ShouldUpdate())
                    {
                        Console.Write("\r{0}, {1:00.0}% complete                        ", p.Status, p.Percent);
                    }
                });
                var run = new PbfLcMsRun(rawFilePath, reader, pbfFilePath, 0, 0, progress, false, startScan, endScan);
                Console.WriteLine();
            }

            Console.WriteLine("PbfFormatVersion: {0}", PbfLcMsRun.FileFormatId);
            return 0;

#if (!DEBUG)
            }
            catch (Exception ex)
            {
                // NOTE: The DMS Analysis Manager looks for this text; do not change it
                Console.WriteLine("Exception while processing: " + ex.Message);
                Console.WriteLine(ex.StackTrace);
                var errorCode = -Math.Abs(ex.Message.GetHashCode());
                if (errorCode == 0)
                    return -1;
                else
                    return errorCode;
            }
#endif
        }

        private static void PrintUsageInfo(string errorMessage = null)
        {
            if (!string.IsNullOrWhiteSpace(errorMessage))
            {
                Console.WriteLine("----------------------------------------------------------");
                Console.WriteLine("Error: " + errorMessage);
                Console.WriteLine("----------------------------------------------------------");
                Console.WriteLine();
            }

            var exePath = System.Reflection.Assembly.GetEntryAssembly().Location;
            var exeName = Path.GetFileName(exePath);

            Console.WriteLine(
                Path.GetFileNameWithoutExtension(exeName) + " " + Version + "\n" +
                "Usage: " + exeName + "\n" +
                "\t-s RawFilePath (*.raw or directory)\n" +
                "\t[-o OutputDir]\n" +
                "\t[-start #]\n" +
                "\t[-end #]\n"
                );

            Console.WriteLine("Optionally use -start and -end to limit the scan range to include in the .pbf file");
            Console.WriteLine("For example " + exeName + " Dataset.raw -start 2000 -end 3000");

            // Wait for 1.5 seconds
            System.Threading.Thread.Sleep(1500);
        }
    }
}
