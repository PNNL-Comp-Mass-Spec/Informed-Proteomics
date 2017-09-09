using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using PRISM;

namespace MSPathFinderT
{
    public class Program
    {
        public const string Name = "MSPathFinderT";
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
            var errorCode = 0;

#if (!DEBUG)
            try
#endif
            {
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, EnableExtendedFlags);

                var entryAsmName = System.Reflection.Assembly.GetEntryAssembly().GetName().Name;
                var parser = new CommandLineParser<TopDownInputParameters>(entryAsmName, Version);
                parser.UsageExamples.Add(string.Format("Perform a target+decoy search with default parameters and a mods file:\n{0}.exe -s testfile.raw -d sampleProteins.fasta -mod mods.txt -tda 1", entryAsmName));
                // TODO: add more examples, maybe a link to GitHub?
                var results = parser.ParseArgs(args);
                if (!results.Success)
                {
                    // Wait for 1.5 seconds
                    System.Threading.Thread.Sleep(1500);
                    return -3;
                }

                if (!results.ParsedResults.Validate())
                {
                    parser.PrintHelp();

                    // Wait for 1.5 seconds
                    System.Threading.Thread.Sleep(1500);

                    return -3;
                }
                var parameters = results.ParsedResults;

                Console.WriteLine(Name + " " + Version);
                parameters.Display();

                parameters.MaxNumNTermCleavages = 1; // max number of N-term cleavages
                parameters.MaxNumCTermCleavages = 0; // max number of C-term cleavages

                foreach (var specFilePath in parameters.SpecFilePaths)
                {
                    // Update the spec file path in the parameters for each search
                    parameters.SpecFilePath = specFilePath;
                    parameters.Write();

                    var topDownLauncher = new IcTopDownLauncher(parameters);

                    var success = topDownLauncher.RunSearch();

                    if (success)
                    {
                        continue;
                    }

                    // topDownLauncher returned false (not successful)

                    // NOTE: The DMS Analysis Manager looks for this text; do not change it
                    var errorMsg = "Error processing " + Path.GetFileName(specFilePath) + ": ";

                    if (string.IsNullOrWhiteSpace(topDownLauncher.ErrorMessage))
                    {
                        errorMsg += "unknown error";
                    }
                    else
                    {
                        errorMsg += topDownLauncher.ErrorMessage;
                    }

                    Console.WriteLine(errorMsg);

                    if (errorCode == 0)
                    {
                        // This is the first error encountered; update the error code
                        // (though we will continue processing the next file if there is one)
                        errorCode = -Math.Abs(errorMsg.GetHashCode());
                        if (errorCode == 0)
                            return -1;
                        else
                            return errorCode;
                    }
                }
            }
#if (!DEBUG)
            catch (Exception ex)
            {
                // NOTE: The DMS Analysis Manager looks for this text; do not change it
                Console.WriteLine("Exception while processing: " + ex.Message);
                Console.WriteLine(ex.StackTrace);
                errorCode = -Math.Abs(ex.Message.GetHashCode());

                System.Threading.Thread.Sleep(1500);

                if (errorCode == 0)
                    return -1;
                else
                    return errorCode;
            }
#endif

                return errorCode;
        }
    }
}
