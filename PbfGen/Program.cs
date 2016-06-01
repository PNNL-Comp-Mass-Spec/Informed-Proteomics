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
        public const string Name = "PbfGen";
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

            try
            {
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, EnableExtendedFlags);

                if (args.Length == 0)
                {
                    PrintUsageInfo();
                    return -1;
                }

                var paramDic = new Dictionary<string, string>
                {
                    {"-s", null},
                    {"-o", null}
                };

                for (var i = 0; i < args.Length / 2; i++)
                {
                    var key = args[2 * i];
                    var value = args[2 * i + 1];
                    if (!paramDic.ContainsKey(key))
                    {
                        PrintUsageInfo("Invalid parameter: " + key);
                        return -1;
                    }
                    paramDic[key] = value;
                }

                // Parse command line parameters
                specFilePath = paramDic["-s"];

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
                outputDir = paramDic["-o"] ?? (specPathIsDirectory ? specFilePath : Path.GetDirectoryName(Path.GetFullPath(specFilePath)));
                if (outputDir == null)
                {
                    PrintUsageInfo("Invalid output file directory: " + specFilePath);
                    return -1;
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
                string[] specFilePaths = new[] { specFilePath };
                if (Directory.Exists(specFilePath) && !MassSpecDataReaderFactory.IsADirectoryDataset(specFilePath))
                {
                    specFilePaths = Directory.GetFiles(specFilePath, "*.raw");
                }

                foreach (var rawFilePath in specFilePaths)
                {
                    var pbfFilePath = Path.Combine(outputDir, Path.GetFileNameWithoutExtension(rawFilePath) + PbfLcMsRun.FileExtension);

                    bool isCurrent;
                    if (File.Exists(pbfFilePath) && PbfLcMsRun.CheckFileFormatVersion(pbfFilePath, out isCurrent) && isCurrent)
                    {
                        Console.WriteLine("{0} already exists.", pbfFilePath);
                        continue;
                    }

                    Console.WriteLine("Creating {0} from {1}", pbfFilePath, rawFilePath);
                    IMassSpecDataReader reader = MassSpecDataReaderFactory.GetMassSpecDataReader(rawFilePath);
                    var progress = new Progress<ProgressData>(p =>
                    {
                        p.UpdateFrequencySeconds = 2;
                        if ((p.Percent % 25).Equals(0) || p.ShouldUpdate())
                        {
                            Console.Write("\r{0}, {1:00.0}% complete                        ", p.Status, p.Percent);
                        }
                    });
                    var run = new PbfLcMsRun(rawFilePath, reader, pbfFilePath, 0, 0, progress);
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

            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t-s RawFilePath (*.raw or directory)\n" +
                "\t[-o OutputDir]\n"
                );

            // Wait for 1.5 seconds
            System.Threading.Thread.Sleep(1500);
        }
    }
}
