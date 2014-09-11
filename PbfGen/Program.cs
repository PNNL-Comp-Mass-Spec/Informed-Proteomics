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
                {"-o", null},
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
            if (!File.Exists(specFilePath))
            {
                PrintUsageInfo("File not found: " + specFilePath + ".");
                return;
            }
            if (Path.GetExtension(specFilePath).ToLower().Equals("raw"))
            {
                PrintUsageInfo("Invalid file extension: " + specFilePath + ".");
                return;
            }
            var outputDir = paramDic["-o"] ?? Path.GetDirectoryName(specFilePath);
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

            var rafFilePath = outputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(specFilePath) + PbfLcMsRun.FileExtension;

            Console.WriteLine("Creating {0} from {1}", rafFilePath, specFilePath);
            var reader = new XCaliburReader(specFilePath);
            var run = new LcMsRun(reader, 0, 0);
            run.WriteAsPbf(rafFilePath);
            Console.WriteLine("RafFormatVersion: {0}", PbfLcMsRun.FileFormatId);
        }

        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine(
                "Usage: " + Name + ".exe\n" +
                "\t-s RawFilePath (*.raw)\n" +
                "\t[-o OutputDir]\n"
                );
        }
    }
}
