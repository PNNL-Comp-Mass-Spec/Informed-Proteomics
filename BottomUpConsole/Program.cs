using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using InformedProteomics.BottomUp.Execution;

namespace MSPathFinder
{
    internal class Program
    {
        public const string Name = "MSPathFinder";
        public const string Version = "0.15 (Sep 11, 2014)";
        public const double CorrThreshold = 0.7;
        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);

        private const uint EnableExtendedFlags = 0x0080;

        private static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length % 2 != 0)
            {
                PrintUsageInfo("The number of arguments must be even.");
                return;
            }

            // initialize parameters
            var paramDic = new Dictionary<string, string>
            {
                {"-s", null},
                {"-d", null},
                {"-o", null},
                {"-e", "1"},
                {"-ntt", "2"},
                {"-mod", null},
                {"-t", "10"},
                {"-f", "10"},
                {"-tda", "0"},
                {"-minLength", "6"},
                {"-maxLength", "40"},
                {"-minCharge", "1"},
                {"-maxCharge", "4"},
                {"-minFragCharge", "1"},
                {"-maxFragCharge", "3"}
            };

            for (var i = 0; i < args.Length / 2; i++)
            {
                var key = args[2 * i];
                var value = args[2 * i + 1];
                if (!paramDic.ContainsKey(key))
                {
                    PrintUsageInfo("Invalid parameter: " + key);
                    return;
                }
                paramDic[key] = value;
            }

            var parameters = new BottomUpInputParameters();
            var message = parameters.Parse(paramDic);
            if (message != null)
            {
                PrintUsageInfo(message);
                return;
            }

            Console.WriteLine(Name + " " + Version);
            parameters.Display();
            parameters.Write();

            foreach (string specFilePath in parameters.SpecFilePaths)
            {
                var bottomUpLauncher = new IcBottomUpLauncher(
                    specFilePath,
                    parameters.DatabaseFilePath,
                    parameters.OutputDir,
                    parameters.AminoAcidSet,
                    parameters.Enzyme,
                    parameters.MinSequenceLength,
                    parameters.MaxSequenceLength,
                    parameters.MinPrecursorIonCharge,
                    parameters.MaxPrecursorIonCharge,
                    parameters.MinProductIonCharge,
                    parameters.MaxProductIonCharge,
                    parameters.PrecursorIonTolerancePpm,
                    parameters.ProductIonTolerancePpm,
                    parameters.Tda,
                    parameters.NumTolerableTermini
                    );

                bottomUpLauncher.RunSearch(CorrThreshold);
            }
        }


        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t-s SpectrumFile or Folder (*.raw)\n" +
                "\t-d DatabaseFile (*.fasta or *.fa)\n" +
                "\t[-o OutputFolder]\n" +
                "\t[-e EnzymeId] (0: Unspecific cleavage, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: AlphaLP, 9: No cleavage)\n" +
                "\t[-ntt NumTolerableTermini] (0, 1, or 2, Default: 2)\n" +
                "\t[-mod ModificationFileName] (modification file, default: no modification)\n" +
                "\t[-t PrecursorToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-f FragmentIonToleranceInPpm] (e.g. 10, Default: 10)\n" +
                "\t[-tda 0/1] (0: don't search decoy database (default), 1: search shuffled decoy database)\n" +
                "\t[-minLength MinSequenceLength] (minimum sequence length, default: 6)\n" +
                "\t[-maxLength MaxSequenceLength] (maximum sequence length, default: 40)\n" +
                "\t[-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 1)\n" +
                "\t[-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 4)\n" +
                "\t[-minFragCharge MinPrecursorCharge] (minimum fragment ion charge, default: 1)\n" +
                "\t[-maxFragCharge MaxPrecursorCharge] (maximum fragment ion charge, default: 3)\n"
                );
        }
    }
}
