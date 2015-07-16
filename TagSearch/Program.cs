using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
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
            const string dataSetPath = @"E:\Jungkap\CompRef";
            const string fastaFilePath = @"E:\Jungkap\CompRef\ID_003278_4B4B3CB1.fasta";
            const string modsFilePath = @"E:\Jungkap\CompRef\Mods.txt";

            var fileEntries = Directory.GetFiles(dataSetPath);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("pbf") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();
            
            var fastaDb = new FastaDatabase(fastaFilePath);
            var tolerance = new Tolerance(10);
            var aaSet = new AminoAcidSet(modsFilePath);

            var dataIndex = new int[3][];
            dataIndex[0] = new int[] {0, 1, 2};
            dataIndex[1] = new int[] { 3, 4, 5 };
            dataIndex[2] = new int[] { 6, 7, 8 };

            //for (var i = 0; i < dataset.Count; i++)
            var setIndex = int.Parse(_paramDic["-i"]);
            foreach(var i in dataIndex[setIndex - 1])
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", dataSetPath, dataset[i]);
                //var ms1File = string.Format(@"{0}\{1}.ms1ft", dataSetPath, dataset[i]);
                var tagFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, ".seqtag");
                var retFile = string.Format(@"{0}\{1}_tagmatch.tsv", dataSetPath, dataset[i]);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                const int minTagLength = 5;
                var tagParser = new SequenceTagParser(tagFilePath, minTagLength, 100);

                Console.WriteLine("Start processing: {0}", rawFile);

                //TestTagBasedSearch(run, tagParser, fastaDb, tolerance, aaSet);
                var engine = new ScanBasedTagSearchEngine(run, tagParser, fastaDb, tolerance, aaSet);
                engine.outWriter = new StreamWriter(retFile);
                engine.RunSearch();

                engine.outWriter.Close();
                
                Console.WriteLine("Complete processing: {0}", rawFile);
            }            
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
