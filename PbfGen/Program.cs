using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.MassSpecData;

namespace PbfGen
{
    public class Program
    {
        public const string Name = "PbfGen";
        public const string Version = "0.62 (Dec 11, 2014)";

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        public static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length == 0)
            {
                PrintUsageInfo();
                return;
            }

            var paramDic = new Dictionary<string, string>
            {
                {"-s", null},
                {"-o", null}
            };

            for (var i = 0; i < args.Length/2; i++)
            {
                var key = args[2*i];
                var value = args[2*i + 1];
                if (!paramDic.ContainsKey(key))
                {
                    PrintUsageInfo("Invalid parameter: " + key);
                    return;
                }
                paramDic[key] = value;
            }

            // Parse command line parameters
            var specFilePath = paramDic["-s"];

            if (specFilePath == null)
            {
                PrintUsageInfo("Missing required parameter -s!");
                return;
            }

            if (!File.Exists(specFilePath) && !Directory.Exists(specFilePath))
            {
                PrintUsageInfo("File not found: " + specFilePath + ".");
                return;
            }

            if (!Directory.Exists(specFilePath) && !Path.GetExtension(specFilePath).ToLower().Equals(".raw") && !Path.GetExtension(specFilePath).ToLower().Equals(".mzml") && !specFilePath.ToLower().EndsWith(".mzml.gz"))
            {
                PrintUsageInfo("Invalid file extension: (" + Path.GetExtension(specFilePath) + ") " + specFilePath + ".");
                return;
            }

            var outputDir = paramDic["-o"] ?? (Directory.Exists(specFilePath) ? specFilePath : Path.GetDirectoryName(specFilePath));
            if (outputDir == null)
            {
                PrintUsageInfo("Invalid raw file directory: " + specFilePath + ".");
                return;
            }

            if (outputDir[outputDir.Length - 1] == Path.DirectorySeparatorChar) outputDir = outputDir.Remove(outputDir.Length - 1);
            if (!Directory.Exists(outputDir))
            {
                if (File.Exists(outputDir) && !File.GetAttributes(outputDir).HasFlag(FileAttributes.Directory))
                {
                    PrintUsageInfo("OutputDir " + outputDir + " is not a directory!");
                    return;
                }
                Directory.CreateDirectory(outputDir);
            }

            var specFilePaths = Directory.Exists(specFilePath) ? Directory.GetFiles(specFilePath, "*.raw") : new[] { specFilePath };

            foreach (var rawFilePath in specFilePaths)
            {
                var rafFilePath = outputDir + Path.DirectorySeparatorChar +
                                           Path.GetFileNameWithoutExtension(rawFilePath) + PbfLcMsRun.FileExtension;

                if(File.Exists(rafFilePath) && PbfLcMsRun.CheckFileFormatVersion(rafFilePath))
                {
                    Console.WriteLine("{0} already exists.", rawFilePath);
                }
                else
                {
                    Console.WriteLine("Creating {0} from {1}", rafFilePath, rawFilePath);
                    IMassSpecDataReader reader = null;
                    if (Path.GetExtension(rawFilePath).ToLower().Equals(".raw"))
                    {
                        reader = new XCaliburReader(rawFilePath);
                    }
                    else if (Path.GetExtension(rawFilePath).ToLower().Equals(".mzml") ||
                             rawFilePath.ToLower().EndsWith(".mzml.gz"))
                    {
                        reader = new MzMLReader(rawFilePath, false, true);
                    }
                    var run = new InMemoryLcMsRun(reader, 0, 0);
                    run.WriteAsPbf(rafFilePath);
                }
            }
            Console.WriteLine("PbfFormatVersion: {0}", PbfLcMsRun.FileFormatId);
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
