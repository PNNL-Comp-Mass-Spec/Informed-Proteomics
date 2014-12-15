using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using LibSVMsharp;
using LibSVMsharp.Extensions;
using LibSVMsharp.Helpers;

namespace MS1FeatureExtractor
{
    class Program
    {
        public const string Name = "MS1FeatureExtractor";
        public const string Version = "0.1 (Nov 24, 2014)";

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;
        
        static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);

            if (args.Length == 0)
            {
                args = new string[]
                {
                    "-i", @"D:\MassSpecFiles\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw", "-minMass", "3000",
                    "-maxMass", "3100", "-svm" , @"D:\Test\ms1fit_svm_model.txt"
                };
            }

            if (args.Length % 2 != 0)
            {
                PrintUsageInfo("The number of arguments must be even.");
                return;
            }



            // initialize parameters
            var paramDic = new Dictionary<string, string>
            {
                {"-i", null},
                {"-svm", null},
                {"-minCharge", "2"},
                {"-maxCharge", "60"},
                {"-minMass", "3000.0"},
                {"-maxMass", "50000.0"},
                {"-post", "n"},
                {"-tmpTsv", null},
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


            if (paramDic["-post"].Equals("y") || paramDic["-post"].Equals("Y"))
            {
                var tsvFile = paramDic["-tmpTsv"];
                if (tsvFile == null) return;
                var j = tsvFile.LastIndexOf('.');
                var outPath = tsvFile.Substring(0, j);
                Ms1FeaturePostProcessing(tsvFile, outPath, 0.01);
                return;
            }
            


            _minSearchMass = double.Parse(paramDic["-minMass"]);
            _maxSearchMass = double.Parse(paramDic["-maxMass"]);
            _minSearchCharge = Int32.Parse(paramDic["-minCharge"]);
            _maxSearchCharge = Int32.Parse(paramDic["-maxCharge"]);
            
            _inputPath = paramDic["-i"];
            _svmPath = paramDic["-svm"];

            var model = SVM.LoadModel(_svmPath);
            if (model == null)
            {
                Console.WriteLine("Failed to load a SVM model in {0}", _svmPath);
                return;
            }

            _predictor = new Ms1FeaturePredictor(model);
            
            var attr = File.GetAttributes(_inputPath);
            if ((attr & FileAttributes.Directory) == FileAttributes.Directory)
                ProcessDirectory(_inputPath);
            else
            {
                ProcessFile(_inputPath);
            }
        }

        public static void ProcessDirectory(string targetDirectory)
        {
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (var fileName in fileEntries)
                ProcessFile(fileName);

            var subdirectoryEntries = Directory.GetDirectories(targetDirectory);
            foreach (var subdirectory in subdirectoryEntries)
                ProcessDirectory(subdirectory);
        }

        public static void ProcessFile(string path)
        {
            if (path.EndsWith(".raw"))
            {
                var j = path.LastIndexOf('.');
                var outPath = path.Substring(0, j);
                var tmpFile = Ms1FeatureExtract(path, outPath);

                Ms1FeaturePostProcessing(tmpFile, outPath);
            }
        }

        private static List<ChargeLcScanCluster>[] _massBinToClusterMap;
        private static ChargeLcScanMatrix _csm;
        private static double _minSearchMass;
        private static double _maxSearchMass;
        private static int _minSearchCharge;
        private static int _maxSearchCharge;
        private static string _inputPath;
        private static string _svmPath;
        private static int _minSearchMassBin;
        private static int _maxSearchMassBin;
        private static Ms1FeaturePredictor _predictor;

        private const int MaxThreads = 4;

        private static string Ms1FeaturePostProcessing(string tsvFile, string outPath, double probTh = 0.1)
        {
            var tsvParser = new TsvFileParser(tsvFile);
            var tsvData = tsvParser.GetAllData();
            var newTsvFile = string.Format("{0}_ms1ft_{1:0.0000}.tsv", outPath, probTh);

            Console.WriteLine("Start Post-Processing using SVM classfier");

            var featureHeaders = new string[]
            {
                "min_charge", "max_charge", 
                "envelope_corr", "summed_envelope_corr", 
                "hypergeometric_score", "summed_hypergeometric_score",
                "ranksum_score","summed_ranksum_score"
            };

            var outFeatureHeaders = new string[]
            {
                "min_scan_num", "max_scan_num", "min_charge", "max_charge",
                "monoisotopic_mw", "rep_scan_num", "rep_charge", "rep_mz",
                "abundance", "isotopic_envelope", "envelope_corr", "summed_envelope_corr", "probability"
            };

            var tsvWriter = new System.IO.StreamWriter(newTsvFile);

            for (var j = 0; j < outFeatureHeaders.Length; j++)
            {
                tsvWriter.Write(outFeatureHeaders[j]);
                tsvWriter.Write(j < outFeatureHeaders.Length - 1 ? "\t" : "\n");
            }
            int count = 0;
            for (var i = 0; i < tsvParser.NumData; i++)
            {
                var nodes = new SVMNode[featureHeaders.Length];
                for (var j = 0; j < featureHeaders.Length; j++)
                {
                    var value = Double.Parse(tsvData[featureHeaders[j]].ElementAt(i));
                    nodes[j] = new SVMNode(j + 1, value);
                }
                var prob = _predictor.PredictProbability(nodes);
                if (prob > probTh)
                {
                    for (var j = 0; j < outFeatureHeaders.Length - 1; j++)
                    {
                        tsvWriter.Write(tsvData[outFeatureHeaders[j]].ElementAt(i));
                        tsvWriter.Write("\t");
                    }
                    tsvWriter.WriteLine(prob);
                    count++;
                }
            }

            tsvWriter.Close();

            Console.WriteLine("Complete Post-Processing. {0} features have been collected", count);

            return newTsvFile;
        }

        private static string Ms1FeatureExtract(string rawFile, string outFile)
        {
            var tsvFilePath = string.Format("{0}_ms1ft.tmp.tsv", outFile);

            Console.WriteLine("Input File : {0}", rawFile);
            Console.WriteLine("Output tsv File : {0}", tsvFilePath);
            
            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine("Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0.0);


            ChargeLcScanMatrix.Initilize(run, _minSearchCharge, _maxSearchCharge);
            Console.WriteLine("Completed loading MS1 data; Elapsed Time = {0:0.00} [s]", (stopwatch.ElapsedMilliseconds)/1000.0d);

            var comparer = ChargeLcScanMatrix.GetMzComparerWithBinning();
            _minSearchMassBin = comparer.GetBinNumber(_minSearchMass);
            _maxSearchMassBin = comparer.GetBinNumber(_maxSearchMass);

            /*
            var binNumberList = Enumerable.Range(_minSearchMassBin, _maxSearchMassBin - _minSearchMassBin + 1).ToList();
            var rnd = new Random();
            binNumberList = binNumberList.OrderBy(item => rnd.Next()).ToList();* 
            Ms1FeatureExtractThread.SetOutputFile(outFile);
            
            Parallel.ForEach(binNumberList, binNum =>
            {
                var csm = ChargeLcScanMatrix.GetInstance();
                var clusters = csm.GetProbableChargeScanClusters(binNum);
                Ms1FeatureExtractThread.WriteToOutputFile(clusters);
            });

            Console.WriteLine(stopwatch.ElapsedMilliseconds/60000.0d);
            */
            
            /*
            var a = new IntRange[MaxThreads];
            var binsPerThread = (int)Math.Round((double)binNumberList.Count / (double)MaxThreads);
            var threads = new Thread[MaxThreads];
            for (var threadNo = 0; threadNo < MaxThreads; threadNo++)
            {
                var stIdx = 0;
                if (threadNo > 0) stIdx = a[threadNo - 1].Max + 1;
                var edIdx = stIdx + binsPerThread;
                if (threadNo == MaxThreads - 1) edIdx = binNumberList.Count - 1;
                a[threadNo] = new IntRange(stIdx, edIdx);

                var binNums = new List<int>(binNumberList.GetRange(stIdx, edIdx - stIdx + 1));
                var ms1Job = new Ms1FeatureExtractThread(binNums, ChargeLcScanMatrix.GetInstance(), threadNo);
                threads[threadNo] = new Thread(new ThreadStart(ms1Job.Run));
                threads[threadNo].Start();

                while (!threads[threadNo].IsAlive);
                Thread.Sleep(1);
            }

            for (var threadNo = 0; threadNo < MaxThreads; threadNo++) threads[threadNo].Join();
            Console.WriteLine();
            Console.WriteLine("==== Completed MS1 Feater Extracting ====");
            */
            
			_massBinToClusterMap = new List<ChargeLcScanCluster>[_maxSearchMassBin - _minSearchMassBin + 1];
            for (var i = 0; i < _massBinToClusterMap.Length; i++) _massBinToClusterMap[i] = new List<ChargeLcScanCluster>();

            var tsvWriter = new System.IO.StreamWriter(tsvFilePath);
            tsvWriter.WriteLine(ChargeLcScanCluster.GetHeaderString());

            stopwatch.Restart();
            Console.WriteLine("Start extracting MS1 features");
            var totalMassBins = _maxSearchMassBin - _minSearchMassBin + 1;

            var binNumCursor = _minSearchMassBin;

            _csm = ChargeLcScanMatrix.GetInstance();

            for (var binNum = _minSearchMassBin; binNum <= _maxSearchMassBin; binNum++)
            {
                ExtractFeatures(binNum);
                /*
                var clusters = csm.GetProbableChargeScanClusters(binNum);
                foreach (var cluster in clusters)
                {
                    tsvWriter.WriteLine(cluster.ToString());
                }*/

                if (binNum > _minSearchMassBin && (binNum - _minSearchMassBin) % 1000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 60000.0d;
                    var processedBins = binNum - _minSearchMassBin;
                    var remaining = (totalMassBins - processedBins) * (elapsed / processedBins);
                    Console.WriteLine("Processed {0}/{1} mass bins; Elapsed Time = {2:0.0} [min]; Remaining Time = {3:0.0} [min]", processedBins, totalMassBins, elapsed, remaining);
                    FlushOutput(binNumCursor, binNum - 500, tsvWriter);
                    binNumCursor = binNum - 500 + 1;
                }
            }

            FlushOutput(binNumCursor, _maxSearchMassBin, tsvWriter);
            tsvWriter.Close();

            return tsvFilePath;
        }
        
        private static void FlushOutput(int startBin, int endBin, StreamWriter tsvWriter)
        {
            for (var binNum = startBin; binNum <= endBin; binNum++)
            {
                var clusters = _massBinToClusterMap[binNum - _minSearchMassBin];
                foreach (var cluster in clusters.Where(x => x.Active == true))
                {
                    tsvWriter.WriteLine(cluster.ToString());
                }
                // release memory
                _massBinToClusterMap[binNum - _minSearchMassBin].Clear();
            }
        }
        
        private static void ExtractFeatures(int binNum)
        {
            var comparer = ChargeLcScanMatrix.GetMzComparerWithBinning();
            var monoIsotopicMass = comparer.GetMzAverage(binNum);

            var neighborMassBins = new List<int>();
            if (binNum > _minSearchMassBin) neighborMassBins.Add(binNum - 1);
            if (binNum < _maxSearchMassBin) neighborMassBins.Add(binNum + 1);
            /*
            for (var i = -2; i <= 2; i++)
            {
                if (i == 0) continue;
                var neighborIsotopebinNum = comparer.GetBinNumber(monoIsotopicMass + i);
                if (neighborIsotopebinNum >= _minSearchMassBin && neighborIsotopebinNum <= _maxSearchMassBin) neighborMassBins.Add(neighborIsotopebinNum);
            }
            */

            foreach (var cluster in _csm.GetProbableChargeScanClusters(monoIsotopicMass))
            {
                // Before adding it, Check if there are existing neighbor that can be merged
                // search mass range = +-2 [Da]
                var foundNeighbor = false;

                foreach (var neighborBin in neighborMassBins)
                {
                    foreach (var neighborCluster in _massBinToClusterMap[neighborBin - _minSearchMassBin])
                    {
                        var massDiff = Math.Abs(neighborCluster.RepresentativeMass - cluster.RepresentativeMass);
                        //if (neighborCluster.Overlaps(cluster) && (massDiff < 1e-9 || Math.Abs(massDiff - 1) < 1e-9 || Math.Abs(massDiff - 2) < 1e-9))
                        if (neighborCluster.Overlaps(cluster) && (massDiff < 1e-9))
                        {
                            if (neighborCluster.Active)
                            {
                                if (neighborCluster.GetScore(ChargeLcScanScore.EnvelopeCorrelation) + neighborCluster.GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) > cluster.GetScore(ChargeLcScanScore.EnvelopeCorrelation) + cluster.GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed))
                                {
                                    cluster.Active = false;
                                    //neighborCluster.Merge(cluster);
                                }
                                else
                                {
                                    neighborCluster.Active = false;
                                    cluster.Active = true;
                                    //cluster.Merge(neighborCluster);
                                }
                            }
                            else
                            {
                                cluster.Active = false;
                            }

                            foundNeighbor = true;
                            break;
                        }
                    }

                    if (foundNeighbor) break;
                }

                _massBinToClusterMap[binNum - _minSearchMassBin].Add(cluster);
            }
        }
        
        private static void PrintUsageInfo(string message = null)
        {
            if (message != null) Console.WriteLine("Error: " + message);
            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder]\n" +
                "\t[-svm SVM Model File]\n" +
                "\t[-minCharge MinPrecursorCharge] (minimum precursor ion charge, default: 2)\n" +
                "\t[-maxCharge MaxPrecursorCharge] (maximum precursor ion charge, default: 60)\n" +
                "\t[-minMass MinSequenceMassInDa] (minimum sequence mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxSequenceMassInDa] (maximum sequence mass in Da, default: 50000.0)\n"
                );
        }
    }

    internal class Ms1FeatureExtractThread
    {
        private readonly List<int> _massBinNums;
        private readonly ChargeLcScanMatrix _csm;
        private readonly int _threadId;
        //private readonly string _outFile;

        private static string _tsvFilePath;
        private static object _sbSync;

        internal Ms1FeatureExtractThread(List<int> massBinNums, ChargeLcScanMatrix csm, int threadId)
        {
            _massBinNums = massBinNums;
            _csm = csm;
            _threadId = threadId;
            //_outFile = outFile;
        }

        internal static void SetOutputFile(string outFile)
        {
            _tsvFilePath = string.Format("{0}_ms1ft.tsv", outFile);
            _sbSync = new Object();
        }
        
        public static void WriteToOutputFile(IEnumerable<ChargeLcScanCluster> clusters)
        {
            lock (_sbSync) // Enter synchronization block
            {
                var tsvWriter = new StreamWriter(new FileStream(_tsvFilePath, FileMode.Append, FileAccess.Write, FileShare.Read));

                foreach (var cluster in clusters)
                {
                    tsvWriter.WriteLine(cluster.ToString());
                }
                tsvWriter.Close();
            }
        }
        
        // This method that will be called when the thread is started
        public void Run()
        {
            Console.WriteLine("[Thread {0}] Started for Total {1} mass bins", _threadId, _massBinNums.Count);
            var stopwatch = Stopwatch.StartNew();
            var bufClusters = new List<ChargeLcScanCluster>();
            for(var i = 0; i < _massBinNums.Count; i++)
            {
                var clusters = _csm.GetProbableChargeScanClusters(_massBinNums[i]);
                
                bufClusters.AddRange(clusters);

                if (i > 0 && (i % 100 == 0 || i == _massBinNums.Count - 1))
                {
                    WriteToOutputFile(bufClusters);
                    bufClusters.Clear();

                    var elapsed = stopwatch.ElapsedMilliseconds / 60000.0d;
                    var remaining = (_massBinNums.Count - i) * (elapsed / (double)i);
                    Console.WriteLine("[Thread {0}] Processed {1}/{2} mass bins; Elapsed Time = {3:0.0} [min]; Remaining Time = {4:0.0} [min]", _threadId, i, _massBinNums.Count, elapsed, remaining);
                }
            }
        }
    }
}
