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
        public const string Version = "0.62 (Dec 11, 2014)";

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

                if (!File.Exists(specFilePath) && !Directory.Exists(specFilePath))
                {
                    PrintUsageInfo("File not found: " + specFilePath + ".");
                    return -1;
                }

                var types = MassSpecDataReaderFactory.GetMassSpecDataTypeFilterList();
                types.Remove(".pbf");

                if (!Directory.Exists(specFilePath) && !(types.Select(ext => specFilePath.ToLower().EndsWith(ext)).Any()))
                {
                    PrintUsageInfo("Invalid file extension: (" + Path.GetExtension(specFilePath) + ") " + specFilePath + ".");
                    return -1;
                }

                outputDir = paramDic["-o"] ?? (Directory.Exists(specFilePath) ? specFilePath : Path.GetDirectoryName(specFilePath));
                if (outputDir == null)
                {
                    PrintUsageInfo("Invalid raw file directory: " + specFilePath + ".");
                    return -1;
                }

                if (outputDir[outputDir.Length - 1] == Path.DirectorySeparatorChar) outputDir = outputDir.Remove(outputDir.Length - 1);
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
                var specFilePaths = Directory.Exists(specFilePath) ? Directory.GetFiles(specFilePath, "*.raw") : new[] { specFilePath };

                foreach (var rawFilePath in specFilePaths)
                {
                    var pbfFilePath = outputDir + Path.DirectorySeparatorChar +
                                               Path.GetFileNameWithoutExtension(rawFilePath) + PbfLcMsRun.FileExtension;

                    if (File.Exists(pbfFilePath) && PbfLcMsRun.CheckFileFormatVersion(pbfFilePath))
                    {
                        Console.WriteLine("{0} already exists.", rawFilePath);
                    }
                    else
                    {
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
                        //var run = new InMemoryLcMsRun(reader, 0, 0, progress);
                        //Console.WriteLine();
                        //run.WriteAsPbf(rafFilePath, progress);
                        InMemoryLcMsRun.ConvertToPbf(rawFilePath, reader, 0, 0, pbfFilePath, progress);
                        Console.WriteLine();
                    }
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

        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t-s RawFilePath (*.raw or directory)\n" +
                "\t[-o OutputDir]\n"
                );
        }
    }
}
