using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SequenceTag;
using InformedProteomics.Backend.Utils;

namespace SeqTagGen
{
    class Program
    {
        public const string Name = "SeqTagGen";
        public static string Version
        {
            get
            {
                var programVersion = System.Reflection.Assembly.GetExecutingAssembly().GetName().Version;
                return string.Format("version {0}.{1}.{2} (May 07, 2015)", programVersion.Major, programVersion.Minor, programVersion.Build);
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
            if (args.Length % 2 != 0)
            {
                PrintUsageInfo("The number of arguments must be even.");
                return;
            }

            // initialize parameters
            _paramDic = new Dictionary<string, string>
            {
                {"-i", null},
                {"-o", null},
                {"-t", "5"},
                {"-minLen", "5"},
                {"-maxTags", "-1"},
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

            _tolerance = new Tolerance(double.Parse(_paramDic["-t"]));
            _minLen = Int32.Parse(_paramDic["-minLen"]);
            _maxTags = Int32.Parse(_paramDic["-maxTags"]);
            _inputPath = _paramDic["-i"];
            _outFolderPath = _paramDic["-o"];

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

        private static void ProcessDirectory(string targetDirectory)
        {
            // Process the list of files found in the directory. 
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (string fileName in fileEntries)
            {
                if (MsRawFile(fileName))
                {
                    var pbfFilePath = PbfLcMsRun.GetPbfFileName(fileName);
                    if (!File.Exists(pbfFilePath)) ProcessFile(fileName);
                }
                else if (MsPbfFile(fileName)) ProcessFile(fileName);
            }
        }

        private static void ProcessFile(string path)
        {
            var rawFile = path;
            var outFile = MassSpecDataReaderFactory.ChangeExtension(path, "seqtag");
            //var tmpOutFile = Path.ChangeExtension(path, "tmp.seqtag");

            if (_outFolderPath != null)
            {
                if (!Directory.Exists(_outFolderPath)) Directory.CreateDirectory(_outFolderPath);
                outFile = _outFolderPath + @"\" + Path.GetFileName(outFile);
            }

            if (File.Exists(outFile))
            {
                Console.WriteLine("Sequence tag result already exists: {0}", outFile);
                return;
            }

            if (!File.Exists(rawFile))
            {
                Console.WriteLine("Cannot find input file: {0}", rawFile);
                return;
            }

            var tmpOutFile = MassSpecDataReaderFactory.ChangeExtension(outFile, "seqtag.tmp");

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine("Input data : {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var ms2ScanNums = run.GetScanNumbers(2);
            var totalScans = ms2ScanNums.Count;
            
            var tmpWriter = new StreamWriter(tmpOutFile);

            tmpWriter.WriteLine("ScanNum\tSequenceTag\tIsPrefix\tFlankingMass\tRankSumPvalue\tMassRMSE");

            var avgTags = 0;
            var nProcessed = 0;
            foreach (var scanNum in ms2ScanNums)
            {
                var ms2Spec = run.GetSpectrum(scanNum) as ProductSpectrum;
                
                var tagFinder = new SequenceTagFinder(ms2Spec, _tolerance, _minLen);

                var nTags = 0;
                foreach (var tag in tagFinder.FindSequenceTags())
                {
                    var flankingMass = tagFinder.DeconvolutedPeaks[tag[0].Node1].Mass;
                    

                    nTags++;
                    if (tag.Count >= 6 || (nTags < _maxTags && _maxTags > 0))
                    {
                        double[] rmse;
                        var tagStrSet = tag.GetTagStrings(out rmse);
                        for(var t = 0; t < tagStrSet.Length; t++)
                        {
                            tmpWriter.WriteLine("{0}\t{1}\t1\t{2}\t{3}\t{4}", scanNum, tagStrSet[t], flankingMass, tag.Score, rmse[t]);
                            tmpWriter.WriteLine("{0}\t{1}\t0\t{2}\t{3}\t{4}", scanNum, SequenceTag.Reverse(tagStrSet[t]), flankingMass, tag.Score, rmse[t]);
                        }    
                    }
                }

                nProcessed++;
                avgTags += nTags;

                if (nProcessed%50 == 0)
                {
                    //var percentage = (double) nProcessed/(double) totalScans;
                    Console.WriteLine("Processing {0}/{1} spectra; Avg. # of tags = {2:0.00};  Elapsed Time = {3:0.000} sec", nProcessed, totalScans, (double)avgTags / (double)nProcessed, (stopwatch.ElapsedMilliseconds) / 1000.0d);
                }
            }
            tmpWriter.Close();
            File.Move(tmpOutFile, outFile);
            Console.WriteLine("Complete sequence tag generation; Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine("Output : {0}", outFile);
        }
        
        private static bool MsRawFile(string path)
        {
            return (path.EndsWith(".raw") || path.EndsWith(".mzML"));
        }
        
        private static bool MsPbfFile(string path)
        {
            return path.EndsWith(".pbf");
        }

        private static Tolerance _tolerance;
        private static int _minLen;
        private static string _inputPath;
        private static string _outFolderPath;
        private static int _maxTags;
        private static Dictionary<string, string> _paramDic;

        private static void PrintUsageInfo(string errorMessage = null)
        {

            if (!string.IsNullOrWhiteSpace(errorMessage))
            {
                Console.WriteLine(@"----------------------------------------------------------");
                Console.WriteLine(@"Error: " + errorMessage);
                Console.WriteLine(@"----------------------------------------------------------");
                Console.WriteLine();
            }

            Console.WriteLine("****** {0}\t{1} ************", Name, Version);
            Console.WriteLine(
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                "\t[-o OutFolder (default: InputFolder)]\n" +
                "\t[-t Tolerance (default: 5 ppm)]\n"+
                "\t[-minLen MinSequenceTagLength] (minimum length of sequence tag, default: 5)\n" +
                "\t[-maxTags MaxNumberOfSequenceTags] (maximum number of sequence tags per spectrum, default: -1, unlimited: -1)\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"
                );
        }

    }
}
