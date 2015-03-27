using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Dynamic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

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
                return string.Format("version {0}.{1}.{2} (March 23, 2015)", programVersion.Major, programVersion.Minor, programVersion.Build);
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

            //_probabilityThreshold = double.Parse(_paramDic["-minProbability"]);
            _minSearchMass = Math.Max(double.Parse(_paramDic["-minMass"]), 3000);
            _maxSearchMass = Math.Min(double.Parse(_paramDic["-maxMass"]), 50000);
            _minSearchCharge = (int)Math.Max(double.Parse(_paramDic["-minCharge"]), 2);
            _maxSearchCharge = (int)Math.Min(double.Parse(_paramDic["-maxCharge"]), 60);
            _inputPath = _paramDic["-i"];
            _maxThreads = Int32.Parse(_paramDic["-maxThreads"]);
       
            _scoreReport = Str2Bool(_paramDic["-score"]);
            _csvOutput  = Str2Bool(_paramDic["-csv"]);
            _tmpOutput = Str2Bool(_paramDic["-tmp"]);
            Console.WriteLine("************ {0}\t{1} ************", Name, Version);
            foreach (var paramName in _paramDic.Keys)
            {
                if (paramName.Equals("-csv") && !_csvOutput) continue;
                if (paramName.Equals("-tmp") && !_tmpOutput) continue;
                Console.WriteLine("{0}\t{1}", paramName, _paramDic[paramName]);
            }
            var attr = File.GetAttributes(_inputPath);


            if ((attr & FileAttributes.Directory) == FileAttributes.Directory)
            {
                //var allFiles = new List<string>();
                ProcessDirectory(_inputPath);
                //foreach (var path in FilterFiles(allFiles)) ProcessFile(path);
            }
            else
            {
                if (!MsRawFile(_inputPath) && !MsPbfFile(_inputPath))
                {
                    Console.WriteLine("Not supported file extension");
                }
                else
                {
                    ProcessFile(_inputPath);    
                }
                
            }
        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }

        private static void ProcessDirectory(string targetDirectory)
        {
            // Process the list of files found in the directory. 
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (string fileName in fileEntries)
            {
                if (MsRawFile(fileName))
                {
                    var pbfFilePath = Path.ChangeExtension(fileName, "pbf");
                    if (!File.Exists(pbfFilePath)) ProcessFile(fileName);
                }
                else if (MsPbfFile(fileName)) ProcessFile(fileName);
            }

            // Recurse into subdirectories of this directory. 
            /*
            var subdirectoryEntries = Directory.GetDirectories(targetDirectory);
            foreach (string subdirectory in subdirectoryEntries)
                ProcessDirectory(subdirectory);
            */
        }

        private static bool MsRawFile(string path)
        {
            return (path.EndsWith(".raw") || path.EndsWith(".mzML"));
        }
        private static bool MsPbfFile(string path)
        {
            return path.EndsWith(".pbf");
        }
        private static void ProcessFile(string path)
        {
            var rawFile = path;
            var outFile = Ms1FeatureMatrix.GetFeatureFilePath(rawFile);

            if (File.Exists(outFile))
            {
                Console.WriteLine("ProMex output already exists: {0}", outFile);
                return;
            }

            if (!File.Exists(rawFile))
            {
                Console.WriteLine("Cannot find input file: {0}", rawFile);
                return;                    
            }

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine("Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile, path.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun);
            var csm = new Ms1FeatureMatrix(run, _minSearchCharge, _maxSearchCharge, _maxThreads);
            Console.WriteLine("Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds)/1000.0d);
            var outputFile = csm.GenerateFeatureFile(rawFile, _minSearchMass, _maxSearchMass, _scoreReport, _csvOutput, _tmpOutput);
        }

        private static double _minSearchMass;
        private static double _maxSearchMass;
        private static int _minSearchCharge;
        private static int _maxSearchCharge;
        private static string _inputPath;
        private static bool _scoreReport;
        private static bool _csvOutput;
        private static bool _tmpOutput;
        private static Dictionary<string, string> _paramDic;
        private static int _maxThreads;

      
        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine("****** {0}\t{1} ************", Name, Version);
            Console.WriteLine(
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                //"\t[-minProbability 0.1 (default: 0.00)]\n" +
                "\t[-minCharge MinCharge] (minimum charge state, default: 2)\n" +
                "\t[-maxCharge MaxCharge] (maximum charge state, default: 60)\n" +
                "\t[-minMass MinMassInDa] (minimum mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)\n" + 
                "\t[-score n (default: n)]\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"
                );
        }

    }
    
}
