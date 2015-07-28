using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.Backend.Utils;

namespace ProMexAlign
{
    class Program
    {
        public const string Name = "ProMexAlign";
        public static string Version
        {
            get
            {
                var programVersion = System.Reflection.Assembly.GetExecutingAssembly().GetName().Version;
                return string.Format("version {0}.{1}.{2} (" + Misc.GetBuildDateTextFromVersion() + ")", programVersion.Major, programVersion.Minor, programVersion.Build);
            }
        }

        private static Dictionary<string, string> _paramDic;

        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);
        private const uint EnableExtendedFlags = 0x0080;

        static int Main(string[] args)
        {
            try
            {
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, EnableExtendedFlags);

                //args = new string[] {"-i", @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf", "-minMass", "3000", "-maxMass", "30000"};

                if (args.Length == 0)
                {
                    PrintUsageInfo();
                    return -1;
                }

                if (args.Length % 2 != 0)
                {
                    PrintUsageInfo("The number of arguments must be even.");
                    return -1;
                }

                // initialize parameters
                _paramDic = new Dictionary<string, string>
                {
                    {"-i", null},
                    {"-o", null},
                    {"-maxThreads", "0"},
                };

                for (var i = 0; i < args.Length / 2; i++)
                {
                    var key = args[2 * i];
                    var value = args[2 * i + 1];
                    if (!_paramDic.ContainsKey(key))
                    {
                        PrintUsageInfo("Invalid parameter: " + key);
                        return -1;
                    }
                    _paramDic[key] = value;
                }

                // Parse command line parameters
                _inputFilePath = _paramDic["-i"];

                if (_inputFilePath == null)
                {
                    PrintUsageInfo("Missing required parameter -i!");
                    return -1;
                }

                if (!LoadInputFile()) return -1;

                _outFilePath = _paramDic["-o"];
                _maxThreads = int.Parse(_paramDic["-maxThreads"]);
                
                if (_outFilePath == null)
                {
                    var outDirectory = Path.GetDirectoryName(Path.GetFullPath(_inputFilePath));
                    _outFilePath = Path.Combine(outDirectory, "promex_alignment.tsv");
                }
                
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception while parsing the command line parameters: " + ex.Message);
                return -5;
            }

#if (!DEBUG)
            try
            {
#endif
                //var param = new Ms1FeatureFinderInputParameter(_paramDic);
                Console.WriteLine("************ {0} {1} ************", Name, Version);
                //param.Display();
                //var launcher = new Ms1FeatureFinderLauncher(param);
                //launcher.Run();
                Run();
#if (!DEBUG)
            }
            catch (Exception ex)
            {
                // NOTE: The DMS Analysis Manager looks for this text; do not change it
                Console.WriteLine("Exception while processing: " + ex.Message);
                Console.WriteLine(ex.StackTrace);
                var errorCode = -Math.Abs(ex.Message.GetHashCode());
                if (errorCode == 0)
                    return -1;
                else
                    return errorCode;
            }
#endif
            return 0;
        }

        private static string _inputFilePath;
        private static string _outFilePath;
        private static int _maxThreads;
        private static TsvFileParser _inputParser;
        private static string[] _colNames = new string[] { "PdbFilePath", "Ms1ftFilePath", "MsPathFinder", "MsAlign"};

        private static bool LoadInputFile()
        {
            if (!File.Exists(_inputFilePath) && !Directory.Exists(_inputFilePath))
            {
                PrintUsageInfo("File not found: " + _inputFilePath);
                return false;
            }

            _inputParser = new TsvFileParser(_inputFilePath);

            for (var i = 0; i < 2; i++)
            {
                if (_inputParser.GetData(_colNames[i]) == null)
                {
                    PrintUsageInfo(string.Format("Cannot find {0} column in {1} ", _inputFilePath, _colNames[i]));
                    return false;
                }                
            }

            return true;
        }

        private static void Run()
        {
            var rawFiles = _inputParser.GetData(_colNames[0]);
            var ms1FtFiles = _inputParser.GetData(_colNames[1]);
            var mspfFiles = _inputParser.GetData(_colNames[2]);
            var msalignFiles = _inputParser.GetData(_colNames[3]);

            var nDatasets = _inputParser.NumData;
            var prsmContainer = new PrSmContainer[nDatasets];

            var dataset = new string[nDatasets];

            for (var i = 0; i < nDatasets; i++)
            {
                var rawFile = rawFiles[i];
                var ms1File = ms1FtFiles[i];

                var dataName = Path.GetFileName(MassSpecDataReaderFactory.RemoveExtension(rawFile));
                dataset[i] = dataName;

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Error: File not found: {0}", rawFile);
                    return;
                }

                if (!File.Exists(ms1File))
                {
                    Console.WriteLine(@"Error: File not found: {0}", ms1File);
                    return;
                }

                rawFiles.Add(rawFile);
                ms1FtFiles.Add(ms1File);

                prsmContainer[i] = new PrSmContainer(i, dataName);

                // load identification results
                //Console.WriteLine(dataset[i]);

                if (mspfFiles != null)
                {
                    var path = mspfFiles[i];
                    if (!File.Exists(path))
                    {
                        Console.WriteLine(@"Error: File not found: {0}", path);
                        return;
                    }
                    prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsPathFinder);                    
                }

                if (msalignFiles != null)
                {
                    var path = msalignFiles[i];
                    if (!File.Exists(path))
                    {
                        Console.WriteLine(@"Error: File not found: {0}", path);
                        return;
                    }
                    prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);                    
                }
                Console.WriteLine("{0} dataset is loaded", dataName);
                Console.WriteLine("Total # of Protein-Spectrum Mathces = {1}, # of unique proteins = {2}", prsmContainer[i].CountIdentifiedScans(), prsmContainer[i].CountIdentifiedUniqueProteoforms());
            }

            var align = new LcMsFeatureAlignment(ms1FtFiles, rawFiles);
            align.AlignFeatures();
            Console.WriteLine("# of aligned features = {0} ", align.CountAlignedFeatures);
            OutputAlignedFeatures(Path.ChangeExtension(_outFilePath, "tmp.tsv"), align, prsmContainer);

            Console.WriteLine("Start refining Abundnace values");
            align.RefineAbundance();
            Console.WriteLine("Complete refining Abundnace values");
            //var alignedFeatureList = align.GetAlignedFeatures();
            OutputAlignedFeatures(_outFilePath, align, prsmContainer);
            Console.WriteLine("Complete alignment task");

            /*
            var writer = new StreamWriter(_outFilePath);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}", dataset[i]);

            writer.Write("\tMSPathFinder_Protein");
            writer.Write("\tMSPathFinder_Sequence");

            writer.Write("\tMSAlign_Protein");
            writer.Write("\tMSAlign_Sequence");

            writer.Write("\n");
            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = GetMinMaxNet(features);
                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }

                var prsmRet = FindProteinSpectrumMatch(prsmContainer, features);
                for (var k = 0; k < 2; k++)
                {
                    writer.Write("\t");
                    writer.Write((prsmRet[k] != null) ? prsmRet[k].ProteinName : "");
                    writer.Write("\t");
                    writer.Write((prsmRet[k] != null) ? prsmRet[k].SequenceText : "");
                }
                writer.Write("\n");
            }
            writer.Close();            */
        }

        private static void OutputAlignedFeatures(string outFilePath, LcMsFeatureAlignment align, PrSmContainer[] prsmContainer = null)
        {
            var writer = new StreamWriter(outFilePath);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");

            for (var i = 0; i < align.CountDatasets; i++)
            {
                var dataName = Path.GetFileName(MassSpecDataReaderFactory.RemoveExtension(align.RawFileList[i]));
                writer.Write("\t{0}", dataName);
            }

            if (prsmContainer != null)
            {
                writer.Write("\tMSPathFinder_Protein");
                writer.Write("\tMSPathFinder_Sequence");
                writer.Write("\tMSAlign_Protein");
                writer.Write("\tMSAlign_Sequence");                
            }

            writer.Write("\n");

            var alignedFeatureList = align.GetAlignedFeatures();

            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = GetMinMaxNet(features);
                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }

                if (prsmContainer != null)
                {
                    var prsmRet = FindProteinSpectrumMatch(prsmContainer, features);
                    for (var k = 0; k < 2; k++)
                    {
                        writer.Write("\t");
                        writer.Write((prsmRet[k] != null) ? prsmRet[k].ProteinName : "");
                        writer.Write("\t");
                        writer.Write((prsmRet[k] != null) ? prsmRet[k].SequenceText : "");
                    }                    
                }

                writer.Write("\n");
            }
            writer.Close();                     
        }



        private static ProteinSpectrumMatch[] FindProteinSpectrumMatch(PrSmContainer[] prsmContainer, LcMsFeature[] alignedFeatures)
        {
            var ret = new ProteinSpectrumMatch[2];
            var allPrsmList = new List<ProteinSpectrumMatch>();
            for (var i = 0; i < prsmContainer.Length; i++)
            {
                if (alignedFeatures[i] == null) continue;

                var prsmList = prsmContainer[i].FindByFeature(alignedFeatures[i], new Tolerance(10));
                if (prsmList.Count < 1)
                {
                    prsmList = prsmContainer[i].FindByFeature(alignedFeatures[i], new Tolerance(13));
                }

                foreach (var prsm in prsmList)
                {
                    if (!allPrsmList.Contains(prsm)) allPrsmList.Add(prsm);
                }
                //allPrsmList.AddRange(prsmList);
            }

            if (allPrsmList.Count < 1) return ret;

            allPrsmList.Sort();
            foreach (var prsm in allPrsmList)
            {
                if (ret[0] == null && prsm.SearchToolType == ProteinSpectrumMatch.SearchTool.MsPathFinder)
                {
                    ret[0] = prsm;
                }
                else if (ret[1] == null && prsm.SearchToolType == ProteinSpectrumMatch.SearchTool.MsAlign)
                {
                    ret[1] = prsm;
                }

                if (ret[0] != null && ret[1] != null) break;
            }

            return ret;
        }

        private static Tuple<double, int, double, double> GetMinMaxNet(IEnumerable<LcMsFeature> features)
        {
            //var minNet = 1.0d;
            //var maxNet = 0d;
            var minElutionTime = double.MaxValue;
            var maxElutionTime = 0d;
            var massList = new List<double>();

            var charge = 0;
            foreach (var f in features)
            {
                if (f == null) continue;
                //minNet = Math.Min(minNet, f.MinNet);
                //maxNet = Math.Max(maxNet, f.MaxNet);
                minElutionTime = Math.Min(minElutionTime, f.MinElutionTime);
                maxElutionTime = Math.Max(maxElutionTime, f.MaxElutionTime);
                massList.Add(f.Mass);
                charge = f.RepresentativeCharge;
            }
            massList.Sort();
            var mass = massList[(int)(massList.Count * 0.5)];

            //return new Tuple<double, int, double, double>(mass, charge, minNet, maxNet);
            return new Tuple<double, int, double, double>(mass, charge, minElutionTime, maxElutionTime);
        }


        private static void PrintUsageInfo(string message = null)
        {
            if (message != null)
            {
                Console.WriteLine("----------------------------------------------------------");
                Console.WriteLine("Error: " + message);
                Console.WriteLine("----------------------------------------------------------");
                Console.WriteLine();
            }

            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFile]\n" +
                "\t[-o OutputFilePath (default : InputFolder\\promex_alignment.tsv)]\n" +
                "\t[-maxThreads 0 (default: 0 (no limit))]\n"
                );

            // Wait for 1.5 seconds
            System.Threading.Thread.Sleep(1500);
        }
    }
}
