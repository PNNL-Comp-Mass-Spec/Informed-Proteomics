using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Utils;
using InformedProteomics.FeatureFinding;
using PRISM;

namespace ProMex
{
    internal static class Program
    {
        public const string Name = "ProMex";
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

        private static int Main(string[] args)
        {
            try
            {
                // Example text:  ProMex version 1.1.7706 (February 5, 2021)
                // (the build date is computed automatically)

                Console.WriteLine(Name + " " + Version);
                Console.WriteLine();

                var osVersionInfo = new OSVersionInfo();
                if (osVersionInfo.GetOSVersion().IndexOf("windows", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    var handle = Process.GetCurrentProcess().MainWindowHandle;
                    SetConsoleMode(handle, EnableExtendedFlags);
                }

                // Uncomment to debug
                // args = new string[] {"-i", @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf", "-minMass", "3000", "-maxMass", "30000"};

                var parser = new CommandLineParser<ProMexInputParameters>(Name, Version);
                parser.UsageExamples.Add("To create a PNG of the features in an existing ms1ft file " +
                "(requires both a .pbf file and a .ms1ft file):\n\tProMex.exe -i dataset.pbf -ms1ft dataset.ms1ft -featureMap");

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
                Console.WriteLine("Exception while parsing the command line parameters: " + ex.Message);
                System.Threading.Thread.Sleep(1500);
                return -5;
            }
        }

        private static int ProcessFiles(ProMexInputParameters options)
        {
#if (!DISABLE_TRYCATCH)
            try
            {
#endif
                if (options.SourceDatasetPaths.Count == 0)
                {
                    ConsoleMsgUtils.ShowWarning("Source dataset path(s) could not be determined; nothing to do");
                    return -6;
                }

                options.Display();

                var launcher = new LcMsFeatureFinderLauncher(options);
                var errorCount = 0;

                foreach (var datasetFilePath in options.SourceDatasetPaths)
                {
                    options.InputPath = datasetFilePath.FullName;

                    int errorCode;

                    if (string.IsNullOrWhiteSpace(options.ExistingFeaturesFilePath))
                    {
                        errorCode = launcher.Run();
                    }
                    else
                    {
                        errorCode = launcher.CreateFeatureMapImage(options.InputPath, options.ExistingFeaturesFilePath);
                    }

                    if (errorCode != 0)
                        errorCount++;

                    Console.WriteLine();
                }

                return errorCount;

#if (!DISABLE_TRYCATCH)
            }
            catch (Exception ex)
            {
                // NOTE: The DMS Analysis Manager looks for this text; do not change it
                ConsoleMsgUtils.ShowError("Exception while processing", ex);

                var errorCode = -Math.Abs(ex.Message.GetHashCode());
                return errorCode == 0 ? -1 : errorCode;
            }
#endif
        }
    }
}
