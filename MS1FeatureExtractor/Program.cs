using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Runtime.InteropServices;
using InformedProteomics.TopDown.Execution;

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
                return string.Format("version {0}.{1}.{2} (April 01, 2015)", programVersion.Major, programVersion.Minor, programVersion.Build);
            }
        }

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length == 0)
            {
                PrintUsageInfo("Input file(or folder) path must be provided.");
                return;                
            }
            if (args.Length % 2 != 0 )
            {
                PrintUsageInfo("The number of arguments must be even.");
                return;
            }

            // initialize parameters
            _paramDic = new Dictionary<string, string>
            {
                {"-i", null},
                {"-minCharge", "2"},
                {"-maxCharge", "60"},
                {"-minMass", "3000.0"},
                {"-maxMass", "50000.0"},
                {"-score", "n"},
                {"-csv", "n"},
                {"-tmp", "n"},
                //{"-quant", "n"},
                {"-maxThreads", "0"},
            };

            for (var i = 0; i < args.Length / 2; i++)
            {
                var key = args[2 * i];
                var value = args[2 * i + 1];
                if (!_paramDic.ContainsKey(key))
                {
                    PrintUsageInfo("Invalid parameter: " + key);
                    return;
                }
                _paramDic[key] = value;
            }

            var param = new Ms1FeatureFinderInputParameter(_paramDic);
            Console.WriteLine("************ {0}\t{1} ************", Name, Version);
            param.Display();
            var launcher = new Ms1FeatureFinderLauncher(param);
            launcher.Run();
        }

        private static Dictionary<string, string> _paramDic;
        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine("****** {0}\t{1} ************", Name, Version);
            Console.WriteLine(
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                "\t[-minCharge MinCharge] (minimum charge state, default: 2)\n" +
                "\t[-maxCharge MaxCharge] (maximum charge state, default: 60)\n" +
                "\t[-minMass MinMassInDa] (minimum mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)\n" + 
                "\t[-score y or n (default: n)]\n" +
                //"\t[-quant y or n] (quantification purpose, default: n\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"
                );
        }

    }
    
}
