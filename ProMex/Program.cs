using System;
using System.Diagnostics;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Utils;
using InformedProteomics.FeatureFinding;
using PRISM;

namespace ProMex
{
    class Program
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

        [STAThread]
        static int Main(string[] args)
        {
            LcMsFeatureFinderInputParameters parameters;

            try
            {
                var osVersionInfo = new clsOSVersionInfo();
                if (osVersionInfo.GetOSVersion().ToLower().Contains("windows"))
                {
                    var handle = Process.GetCurrentProcess().MainWindowHandle;
                    SetConsoleMode(handle, EnableExtendedFlags);
                }

                //args = new string[] {"-i", @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf", "-minMass", "3000", "-maxMass", "30000"};

                var parser = new CommandLineParser<ProMexInputParameters>(Name, Version);
                parser.UsageExamples.Add("To create a PNG of the features in an existing ms1ft file " +
                "(requires both a .pbf file and a .ms1ft file):\n\tProMex.exe -i dataset.pbf -ms1ft dataset.ms1ft -featureMap");

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

                parameters = results.ParsedResults;
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

                // Example text:  ProMex version 1.0.6527 (November 14, 2017)
                // (the build date is computed automatically)
                Console.WriteLine("************ {0} {1} ************", Name, Version);
                parameters.Display();
                var launcher = new LcMsFeatureFinderLauncher(parameters);
                int errorCode;

                if (string.IsNullOrWhiteSpace(parameters.ExistingFeaturesFilePath))
                {
                    errorCode = launcher.Run();
                }
                else
                {
                    errorCode = launcher.CreateFeatureMapImage(parameters.InputPath, parameters.ExistingFeaturesFilePath);
                }

                return errorCode;
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

                return errorCode;
            }
#endif
        }
    }
}
