using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.Backend.Utils;

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
            {
#endif

            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length == 0)
            {
                PrintUsageInfo();
                return -1;
            }

            if (args.Length % 2 != 0)
            {
                PrintUsageInfo("The number of arguments must be even");
                return -1;
            }

            // initialize parameters
            var paramDic = new Dictionary<string, string>
                {
                    {"-s", null},
                    {"-d", null},
                    {"-o", null},
                    {"-m", "1"},
                    {"-mod", null},
                    {"-t", "10"},
                    {"-f", "10"},
                    {"-tda", "0"},
                    {"-minLength", "21"},
                    {"-maxLength", "300"},
                    {"-minCharge", "2"},
                    {"-maxCharge", "60"},
                    {"-minFragCharge", "1"},
                    {"-maxFragCharge", "20"},
                    {"-minMass", "2000.0"},
                    {"-maxMass", "50000.0"},
                    {"-feature", null},
                    {"-threads", "0"},
                };

            for (var i = 0; i < args.Length / 2; i++)
            {
                var key = args[2 * i];
                var value = args[2 * i + 1];
                if (!paramDic.ContainsKey(key))
                {
                    PrintUsageInfo("Invalid parameter: " + key);
                    return -2;
                }
                paramDic[key] = value;
            }

            var parameters = new TopDownInputParameters();
            var message = parameters.Parse(paramDic);
            if (message != null)
            {
                PrintUsageInfo(message);
                return -3;
            }

            Console.WriteLine(Name + " " + Version);
            parameters.Display();
            parameters.Write();

            foreach (var specFilePath in parameters.SpecFilePaths)
            {
                
                var topDownLauncher = new IcTopDownLauncher(
                    specFilePath,
                    parameters.DatabaseFilePath,
                    parameters.OutputDir,
                    parameters.AminoAcidSet,
                    parameters.MinSequenceLength,
                    parameters.MaxSequenceLength,
                    1, // max number of N-term cleavages
                    0, // max number of C-term cleavages
                    parameters.MinPrecursorIonCharge,
                    parameters.MaxPrecursorIonCharge,
                    parameters.MinProductIonCharge,
                    parameters.MaxProductIonCharge,
                    parameters.MinSequenceMass,
                    parameters.MaxSequenceMass,
                    parameters.PrecursorIonTolerancePpm,
                    parameters.ProductIonTolerancePpm,
                    parameters.Tda,
                    parameters.SearchMode,
                    parameters.FeatureFilePath,
                    parameters.MaxNumThreads
                    );
                
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

#if (!DEBUG)
            }
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
                "\t-s SpectrumFile (*.raw)\n" +
                "\t-d DatabaseFile (*.fasta or *.fa)\n" +
                "\t[-o OutputFolder]\n" +
                "\t[-m SearchMode] (0: multiple internal cleavages, 1: single internal cleavage (default), 2: no internal cleavage)\n" +
                "\t[-mod ModificationFileName] (modification file, default: no modification)\n" +
                "\t[-t PrecursorToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-f FragmentIonToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-tda 0/1] (0: don't search decoy database (default), 1: search shuffled decoy database)\n" +
                "\t[-minLength MinSequenceLength] (minimum sequence length, default: 21)\n" +
                "\t[-maxLength MaxSequenceLength] (maximum sequence length, default: 500)\n" +
                "\t[-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 2)\n" +
                "\t[-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 50)\n" +
                "\t[-minFragCharge MinPrecursorCharge] (minimum fragment ion charge, default: 1)\n" +
                "\t[-maxFragCharge MaxPrecursorCharge] (maximum fragment ion charge, default: 20)\n" +
                "\t[-minMass MinSequenceMassInDa] (minimum sequence mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxSequenceMassInDa] (maximum sequence mass in Da, default: 50000.0)\n" +
                "\t[-feature FeatureFile] (*.ms1ft, *_isos.csv, or *.msalign, default: Run ProMex)\n" +
                "\t[-threads MaxNumThreads] (maximum number of threads, default: 0 automatic setting)\n"
                );

            // Wait for 1.5 seconds
            System.Threading.Thread.Sleep(1500);
        }

    }
}
