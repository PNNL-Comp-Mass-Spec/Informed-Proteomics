using System;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace TopDownConsole
{
    class Program
    {
        public const string Version = "0.1 (June 10, 2014)";
        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            PrintUsageInfo();
            foreach (var s in args)
            {
                Console.WriteLine(s);
            }
        }

        static void PrintUsageInfo()
        {
            Console.WriteLine(
                "POPSICLE " + Version + "\n" +
                "Usage: POPSICLE.exe\n" +
                "\t-s SpectrumFile (*.raw)\n" +
                "\t-d DatabaseFile (*.fasta or *.fa)\n" +
                "\t-o OutputFolder\n" +
                "\t[-m SearchMode] (0: multiple internal cleavages, 1: single internal cleavage (default), 2: no internal cleavage)\n" +
                "\t[-mod ModificationFileName] (modification file, default: no modification)\n" +
                "\t[-t PrecursorToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-f FragmentIonToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-tda 0/1] (0: don't search decoy database (default), 1: search shuffled decoy database)\n" +
                "\t[-minLength MinSequenceLength] (minimum sequence length, default: 21)\n" +
                "\t[-maxLength MaxSequenceLength] (maximum sequence length, default: 300)\n" +
                "\t[-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 2)\n" +
                "\t[-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 30)\n" +
                "\t[-minFragCharge MinPrecursorCharge] (minimum fragment ion charge, default: 1)\n" +
                "\t[-maxFragCharge MaxPrecursorCharge] (maximum fragment ion charge, default: 15)\n" +
                "\t[-minMass MinSequenceMassInDa] (minimum sequence mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxSequenceMassInDa] (maximum sequence mass in Da, default: 50000.0)\n"
                );
        }
    }
}
