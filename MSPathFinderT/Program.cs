using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using PRISM;

namespace MSPathFinderT
{
    public static class Program
    {
        // Ignore Spelling: tda

        public const string Name = "MSPathFinderT";

        /// <summary>
        /// Program version, including the build date
        /// </summary>
        /// <remarks>Example: MSPathFinderT version 1.0.7569 (September 21, 2020)</remarks>
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
            try
            {
                // Example text:  MSPathFinderT version 1.1.7706 (February 5, 2021)
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

                var parser = new CommandLineParser<TopDownInputParameters>(entryAsmName, Version)
                {
                    ParamKeysFieldWidth = 20
                };

                parser.UsageExamples.Add(string.Format("Perform a target+decoy search with default parameters and a mods file:\n{0}.exe -s TestFile.raw -d sampleProteins.fasta -mod mods.txt -tda 1", entryAsmName));

                // TODO: add more examples, maybe a link to GitHub?
                var results = parser.ParseArgs(args);
                if (!results.Success)
                {
                    System.Threading.Thread.Sleep(1500);
                    return -3;
                }

                if (!results.ParsedResults.Validate(results.ParamFilePath))
                {
                    parser.PrintHelp();

                    System.Threading.Thread.Sleep(1500);
                    return -3;
                }

                var options = results.ParsedResults;

                var errorCode = ProcessFiles(results, options);

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

        private static int ProcessFiles(CommandLineParser<TopDownInputParameters>.ParserResults results, TopDownInputParameters options)
        {
#if (!DISABLE_TRYCATCH)
            try
#endif
            {
                options.Display(results.ParamFilePath);

                options.MaxNumNTermCleavages = 1; // max number of N-term cleavages
                options.MaxNumCTermCleavages = 0; // max number of C-term cleavages

                FlipSwitch.Instance.SetUseFlipScoring(options.UseFLIP);

                var errorCode = 0;

                foreach (var specFile in options.SpecFilePaths)
                {
                    // Update the spec file path in the parameters for each search
                    options.SpecFilePath = specFile.FullName;
                    options.Write();

                    var topDownLauncher = new IcTopDownLauncher(options);
                    topDownLauncher.ErrorEvent += TopDownLauncher_ErrorEvent;
                    topDownLauncher.WarningEvent += TopDownLauncher_WarningEvent;
                    topDownLauncher.StatusEvent += TopDownLauncher_StatusEvent;
                    topDownLauncher.ProgressUpdate += TopDownLauncher_ProgressUpdate;

                    var success = topDownLauncher.RunSearch();

                    Console.WriteLine();

                    if (success)
                    {
                        continue;
                    }

                    // topDownLauncher returned false (not successful)

                    // NOTE: The DMS Analysis Manager looks for this text; do not change it
                    var errorMsg = "Error processing " + specFile.Name + ": ";

                    if (string.IsNullOrWhiteSpace(topDownLauncher.ErrorMessage))
                    {
                        errorMsg += "unknown error";
                    }
                    else
                    {
                        errorMsg += topDownLauncher.ErrorMessage;
                    }

                    Console.WriteLine(errorMsg);

                    if (errorCode != 0)
                        continue;

                    errorCode = -Math.Abs(errorMsg.GetHashCode());

                    if (errorCode == 0)
                        errorCode = -1;
                }

                return errorCode;
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

        private static void TopDownLauncher_ProgressUpdate(string progressMessage, float percentComplete)
        {
            Console.WriteLine(progressMessage);
        }

        private static void TopDownLauncher_StatusEvent(string message)
        {
            Console.WriteLine(message);
        }

        private static void TopDownLauncher_ErrorEvent(string message, Exception ex)
        {
            ConsoleMsgUtils.ShowError(message, ex);
        }

        private static void TopDownLauncher_WarningEvent(string message)
        {
            ConsoleMsgUtils.ShowWarning("Warning: " + message);
        }
    }
}
