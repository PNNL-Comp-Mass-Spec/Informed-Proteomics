using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.TagBasedSearch;

namespace TagSearch
{
    class Program
    {
        public const string Name = "TagSearch";
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
                    {"-minMass", "3000.0"},
                    {"-maxMass", "50000.0"},

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
                /*
                var inputFilePath = _paramDic["-i"];
                
                if (inputFilePath == null)
                {
                    PrintUsageInfo("Missing required parameter -i!");
                    return -1;
                }

                if (!File.Exists(inputFilePath) && !Directory.Exists(inputFilePath))
                {
                    PrintUsageInfo("File not found: " + inputFilePath);
                    return -1;
                }
                */
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


        private static void Run()
        {
            const string dataSetPath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            const string fastaFilePath = @"E:\Jungkap\UTEX\db\ID_004352_779A168F.fasta";
            const string modsFilePath = @"E:\Jungkap\UTEX\db\Mods.txt";
            const string seqtagFilePath = @"E:\Jungkap\UTEX";
            var fileEntries = Directory.GetFiles(@"E:\Jungkap\UTEX");

            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();
            
            var fastaDb = new FastaDatabase(fastaFilePath);
            var tolerance = new Tolerance(10);
            var aaSet = new AminoAcidSet(modsFilePath);

            var dataIndex = new int[5][];

            dataIndex[0] = new int[] {0, 1, 2, 3, 4, 5, 6};
            dataIndex[1] = new int[] {7, 8, 9, 10, 11, 12, 13};
            dataIndex[2] = new int[] {14, 15, 16, 17, 18, 19};
            dataIndex[3] = new int[] {20, 21, 22, 23, 24, 25};
            dataIndex[4] = new int[] {26, 27, 28, 29, 30, 31};

            var setIndex = int.Parse(_paramDic["-i"]);
            foreach(var i in dataIndex[setIndex - 1])
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", dataSetPath, dataset[i]);
                var tagFilePath = string.Format(@"{0}\{1}.seqtag", seqtagFilePath, dataset[i]);
                var retFile = string.Format(@"{0}\{1}_tagmatch.tsv", seqtagFilePath, dataset[i]);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                const int minTagLength = 5;
                //var tagParser = new SequenceTagParser(tagFilePath, minTagLength, 100);

                Console.WriteLine("Start processing: {0}", rawFile);
                var tagGen = new SequenceTagGenerator(run, new Tolerance(8));
                var engine = new ScanBasedTagSearchEngine(run, tagGen, new LcMsPeakMatrix(run), fastaDb, tolerance, aaSet);
                
                engine.RunSearch();

                
                
                Console.WriteLine("Complete processing: {0}", rawFile);
            }            
        }

        private static void PrintUsageInfo(string errorMessage = null)
        {
            if (!string.IsNullOrWhiteSpace(errorMessage))
            {
                Console.WriteLine(@"----------------------------------------------------------");
                Console.WriteLine(@"Error: " + errorMessage);
                Console.WriteLine(@"----------------------------------------------------------");
                Console.WriteLine();
            }

            Console.WriteLine(
                Name + " " + Version + "\n" +
                "Usage: " + Name + ".exe\n" +
                "\t[-i InputFolder or InputFile]\n" +
                "\t[-d Fasta file]\n" +
                "\t[-minMass MinMassInDa] (minimum mass in Da, default: 3000.0)\n" +
                "\t[-maxMass MaxMassInDa] (maximum mass in Da, default: 50000.0)\n"
                );

            // Wait for 1.5 seconds
            System.Threading.Thread.Sleep(1500);
        }

    }
    
}
