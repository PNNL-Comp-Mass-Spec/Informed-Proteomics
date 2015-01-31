using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net.NetworkInformation;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using LibSVMsharp;
using LibSVMsharp.Helpers;

namespace Mspot
{
    class Program
    {
        public const string Name = "ProMex";
        public const string Version = "ver. 1.0 (Jan 29, 2014)";

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length == 0)
            {
                PrintUsageInfo();
                return;
            }

            if (args.Length % 2 != 0)
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
                {"-massCollapse", "n"},
                {"-score", "n"},
                {"-csv", "y"},
                {"-minProbability", "0.1"},
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

            _probabilityThreshold = double.Parse(_paramDic["-minProbability"]);
            _minSearchMass = double.Parse(_paramDic["-minMass"]);
            _maxSearchMass = double.Parse(_paramDic["-maxMass"]);
            _minSearchCharge = Int32.Parse(_paramDic["-minCharge"]);
            _maxSearchCharge = Int32.Parse(_paramDic["-maxCharge"]);
            _inputPath = _paramDic["-i"];
            _maxThreads = Int32.Parse(_paramDic["-maxThreads"]);
            
            /*
            var svmFile = paramDic["-svm"];
            var svmModel = SVM.LoadModel(svmFile);
            if (svmModel == null)
            {
                Console.WriteLine("Can not load svm file: {0}", svmFile);
                return;
            }
            _predictor = new Ms1FeatureSvmPredictor(svmModel);
            */

            _scoreReport = Str2Bool(_paramDic["-score"]);
            _massCollapse = Str2Bool(_paramDic["-massCollapse"]);
            _csvOutput  = Str2Bool(_paramDic["-csv"]);
            
            Console.WriteLine("****** {0}\t{1} ************", Name, Version);
            foreach (var paramName in _paramDic.Keys)
                Console.WriteLine("{0}\t{1}", paramName, _paramDic[paramName]);
            
            var attr = File.GetAttributes(_inputPath);
            if ((attr & FileAttributes.Directory) == FileAttributes.Directory)
            {
                var allFiles = new List<string>();
                ProcessDirectory(_inputPath, ref allFiles);
                foreach (var path in FilterFiles(allFiles)) ProcessFile(path);
            }
            else
            {
                ProcessFile(_inputPath);
            }
        }

        private static bool Str2Bool(string value)
        {
            return (value.Equals("y") || value.Equals("Y"));
        }

        public static void ProcessDirectory(string targetDirectory, ref List<string> fileNames)
        {
            //var fileNames = new List<string>();
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (var fileName in fileEntries)
            {
                fileNames.Add(fileName);
            }

            var subdirectoryEntries = Directory.GetDirectories(targetDirectory);
            foreach (var subdirectory in subdirectoryEntries)
                ProcessDirectory(subdirectory, ref fileNames);
        }

        public static void ProcessFile(string path)
        {
            if (path.EndsWith(".raw") || path.EndsWith(".pbf"))
            {
                var rawFile = path;
                var outFile = ChargeLcScanMatrix.GetFeatureFilePath(rawFile);

                if (File.Exists(outFile)) return;

                var stopwatch = Stopwatch.StartNew();
                Console.WriteLine("Start loading MS1 data from {0}", rawFile);
                var run = PbfLcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0.0);
                var csm = new ChargeLcScanMatrix(run, _minSearchCharge, _maxSearchCharge, _maxThreads);
                Console.WriteLine("Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

                var outputFile = csm.GenerateFeatureFile(rawFile, _minSearchMass, _maxSearchMass, _massCollapse, _probabilityThreshold, _scoreReport, _csvOutput);
            }
        }

        private static double _minSearchMass;
        private static double _maxSearchMass;
        private static int _minSearchCharge;
        private static int _maxSearchCharge;
        private static string _inputPath;
        private static bool _massCollapse;
        private static bool _scoreReport;
        private static bool _csvOutput;
        private static double _probabilityThreshold;
        private static IMs1FeaturePredictor _predictor;
        private static Dictionary<string, string> _paramDic;
        private static int _maxThreads;
      
        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                "\t[-minProbability 0.1 (default: 0.1)]\n" +
                "\t[-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 2)\n" +
                "\t[-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 60)\n" +
                "\t[-minMass MinSequenceMassInDa] (minimum sequence mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxSequenceMassInDa] (maximum sequence mass in Da, default: 50000.0)\n" + 
                "\t[-massCollapse n (default: n)]\n" +
                "\t[-score n (default: n)]\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"
                );
        }


        private static List<string> FilterFiles(List<string> fileNames)
        {
            var files = new Dictionary<string, bool>();
            var filteredFiles = new List<string>();

            foreach (var path in fileNames)
            {
                bool pbf = path.EndsWith(".pbf");
                bool raw = path.EndsWith(".raw");
                if (pbf || raw)
                {
                    var j = path.LastIndexOf('.');
                    var outPath = path.Substring(0, j);

                    if (pbf)
                    {
                        if (files.ContainsKey(outPath))
                        {
                            files.Remove(outPath);
                            files.Add(outPath, true);
                        }
                        else
                        {
                            files.Add(outPath, pbf);
                        }
                    }
                    else
                    {
                        if (files.ContainsKey(outPath))
                        {

                        }
                        else
                        {
                            files.Add(outPath, pbf);
                        }
                    }
                }
            }

            foreach (var path in files.Keys)
            {
                filteredFiles.Add(files[path] ? string.Format("{0}.pbf", path) : string.Format("{0}.raw", path));
            }
            return filteredFiles;
        }
    }
    
}
