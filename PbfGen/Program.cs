using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.MassSpecData;
using PRISM;

namespace PbfGen
{
    public class Program
    {
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
            PbfGenInputParameters options;

            try
            {
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, EnableExtendedFlags);

                var exeName = System.Reflection.Assembly.GetEntryAssembly().GetName().Name;

                var parser = new CommandLineParser<PbfGenInputParameters>(exeName, Version);
                parser.UsageExamples.Add($"Using -start and -end to limit the scan range to include in the .pbf file\n\t{exeName}.exe -s Dataset.raw -start 2000 -end 3000");

                var results = parser.ParseArgs(args);

                if (!results.Success)
                {
                    // Wait for 1.5 seconds
                    System.Threading.Thread.Sleep(1500);

                    return -1;
                }

                if (!results.ParsedResults.Validate())
                {
                    parser.PrintHelp();

                    // Wait for 1.5 seconds
                    System.Threading.Thread.Sleep(1500);

                    return -1;
                }

                options = results.ParsedResults;
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception while parsing the command line parameters: " + ex.Message);
                return -5;
            }

#if (!DEBUG)
            try
#endif
            {
                var specFilePaths = new[] {options.SourcePath};
                if (Directory.Exists(options.SourcePath) && !MassSpecDataReaderFactory.IsADirectoryDataset(options.SourcePath))
                {
                    specFilePaths = Directory.GetFiles(options.SourcePath, "*.raw"); // TODO: Support folders with other formats in them too...
                }

                foreach (var rawFilePath in specFilePaths)
                {
                    var pbfFileName = MassSpecDataReaderFactory.ChangeExtension(rawFilePath, PbfLcMsRun.FileExtensionConst);
                    var pbfFilePath = Path.Combine(options.OutputDir, Path.GetFileName(pbfFileName));

                    if (File.Exists(pbfFilePath) && PbfLcMsRun.CheckFileFormatVersion(pbfFilePath, out var isCurrent) && isCurrent)
                    {
                        Console.WriteLine("{0} already exists.", pbfFilePath);
                        continue;
                    }

                    Console.WriteLine("Creating {0} from {1}", pbfFilePath, rawFilePath);

                    if (options.StartScan > 0 && options.EndScan > 0)
                        Console.WriteLine("Only including scans {0} to {1}", options.StartScan, options.EndScan);
                    else if (options.StartScan > 0)
                        Console.WriteLine("Only including scans {0} to the end", options.StartScan);
                    else if (options.EndScan > 0)
                        Console.WriteLine("Only including scans 1 to {0}", options.EndScan);

                    var reader = MassSpecDataReaderFactory.GetMassSpecDataReader(rawFilePath);
                    var progress = new Progress<ProgressData>(p =>
                    {
                        p.UpdateFrequencySeconds = 2;
                        if ((p.Percent % 25).Equals(0) || p.ShouldUpdate())
                        {
                            Console.Write("\r{0}, {1:00.0}% complete                        ", p.Status, p.Percent);
                        }
                    });
                    var run = new PbfLcMsRun(rawFilePath, reader, pbfFilePath, 0, 0, progress, false, options.StartScan, options.EndScan);
                    Console.WriteLine();
                }

                Console.WriteLine("PbfFormatVersion: {0}", PbfLcMsRun.FileFormatId);
                return 0;
            }
#if (!DEBUG)
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
    }
}
