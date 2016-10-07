using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    class TestLcMsFeatureAlignment
    {
        [Test]
        public void TestPeptidomics()
        {
            const string ms1ftFolder = @"\\protoapps\UserData\Jungkap\Mowei\Quant";
            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            string outFilePath = string.Format(@"{0}\aligned_features2.tsv", ms1ftFolder);

            var fileEntries = Directory.GetFiles(ms1ftFolder);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();
            /*var dataset = new String[]
            {
                "MZ20150729FG_WT1",
                "MZ20150729FG_WT2a",
                "MZ20150725FG_WT1",
                "MZ20150725FG_WT2a",
                "MZ20150729FG_MT1",
                "MZ20150729FG_MT2",
                "MZ20150725FG_MT1",
                "MZ20150725FG_MT2"
            };*/

            //var dataset = new List<string>();
            //for (var i = 1; i <= 10; i++) dataset.Add(string.Format(@"CPTAC_Intact_rep{0}_15Jan15_Bane_C2-14-08-02RZ", i));
            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();
            foreach (string datasetName in dataset)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, datasetName);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, datasetName);

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }

                if (!File.Exists(ms1File))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", ms1File);
                    continue;
                }

                rawFiles.Add(rawFile);
                ms1FtFiles.Add(ms1File);
            }

            if (rawFiles.Count == 0)
            {
                Assert.Ignore("No files found.");
            }

            RunFeatureAlignment(ms1FtFiles, rawFiles, outFilePath);
        }

        [Test]
        public void TestCptac10Replicates()
        {
            const string ms1ftFolder = @"D:\MassSpecFiles\CPTAC_rep10";
            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1";
            string outFilePath = string.Format(@"{0}\aligned_features.tsv", ms1ftFolder);

            //var fileEntries = Directory.GetFiles(ms1ftFolder);
            //var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            //dataset.Sort();
            var dataset = new List<string>();
            for (var i = 1; i <= 10; i++) dataset.Add(string.Format(@"CPTAC_Intact_rep{0}_15Jan15_Bane_C2-14-08-02RZ", i));

            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();

            foreach (string datasetName in dataset)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, datasetName);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, datasetName);

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }

                if (!File.Exists(ms1File))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", ms1File);
                    continue;
                }

                rawFiles.Add(rawFile);
                ms1FtFiles.Add(ms1File);
            }

            if (rawFiles.Count == 0)
            {
                Assert.Ignore("No files found.");
            }

            RunFeatureAlignment(ms1FtFiles, rawFiles, outFilePath);
        }

        [Test]
        public void TestCompRef()
        {
            const string outFilePath = @"\\protoapps\UserData\Jungkap\CompRef\aligned\aligned_features_requant.tsv";
            const string ms1ftFolder = @"\\protoapps\UserData\Jungkap\CompRef\ms1ft";
            const string rawFolder = @"\\protoapps\UserData\Jungkap\CompRef\raw";

            var fileEntries = Directory.GetFiles(ms1ftFolder);
            var ms1FtFiles = fileEntries.Where(f => f.EndsWith(".ms1ft")).ToList();

            fileEntries = Directory.GetFiles(rawFolder);
            var rawFiles = fileEntries.Where(f => f.EndsWith(".pbf")).ToList();

            RunFeatureAlignment(ms1FtFiles, rawFiles, outFilePath);
        }

        [Test]
        public void TestIMER()
        {
            const string outFilePath = @"D:\MassSpecFiles\IMER\aligned_features.tsv";
            const string ms1ftFolder = @"D:\MassSpecFiles\IMER";
            const string rawFolder = @"D:\MassSpecFiles\IMER";

            var fileEntries = Directory.GetFiles(ms1ftFolder);
            var ms1FtFiles = fileEntries.Where(f => f.EndsWith(".ms1ft")).ToList();

            fileEntries = Directory.GetFiles(rawFolder);
            var rawFiles = fileEntries.Where(f => f.EndsWith(".pbf")).ToList();

            RunFeatureAlignment(ms1FtFiles, rawFiles, outFilePath);
        }

        [Test]
        public void TestTempCompRefLcMsFeatureAlign()
        {
            const string dataFolder = @"D:\MassSpecFiles\CompRef";
            const string fastaFilePath = @"D:\MassSpecFiles\CompRef\db\ID_003278_4B4B3CB1.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            fastaDb.Read();

            var fileEntries = Directory.GetFiles(dataFolder);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("pbf") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();

            for (var i = 0; i < dataset.Count; i++)
            {
                var writer =
                    new StreamWriter(string.Format(@"D:\MassSpecFiles\CompRef\MsPathFinderMerged\{0}_IcTda.tsv",
                        dataset[i]));

                writer.Write("Scan");
                writer.Write("\t");
                writer.Write("Sequence");
                writer.Write("\t");
                writer.Write("Modifications");
                writer.Write("\t");
                writer.Write("Mass");
                writer.Write("\t");
                writer.Write("ProteinName");
                writer.Write("\t");
                writer.Write("ProteinDesc");
                writer.Write("\t");
                writer.Write("Start");
                writer.Write("\t");
                writer.Write("End");
                writer.Write("\t");
                writer.Write("#MatchedFragments");
                writer.Write("\t");
                writer.Write("QValue");
                writer.Write("\n");

                var path1 = string.Format(@"D:\MassSpecFiles\CompRef\MsPathFinder\{0}_IcTda.tsv", dataset[i]);
                var parser1 = new TsvFileParser(path1);
                OutputMergedResult(writer, parser1, fastaDb);

                var path2 = string.Format(@"D:\MassSpecFiles\CompRef\seqtag\{0}_tagmatch.tsv", dataset[i]);
                var parser2 = new TsvFileParser(path2);
                OutputMergedResult(writer, parser2, fastaDb);
                writer.Close();
            }
        }

        private void OutputMergedResult(StreamWriter writer, TsvFileParser parser, FastaDatabase fastaDb)
        {
            var scoreColumn = parser.GetData("#MatchedFragments") ?? parser.GetData("Score");
            var qValColumn = parser.GetData("QValue");

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Sequence")[i];
                var scanNum = int.Parse(parser.GetData("Scan")[i]);
                var mass = double.Parse(parser.GetData("Mass")[i]);
                var protName = parser.GetData("ProteinName")[i];
                var protDesc = fastaDb.GetProteinDescription(protName);

                var firstResId = int.Parse(parser.GetData("Start")[i]);
                var lastResId = int.Parse(parser.GetData("End")[i]);
                var score = double.Parse(scoreColumn[i]);
                var mod = parser.GetData("Modifications")[i];
                var qvalue = (qValColumn != null) ? qValColumn[i] : "0";

                writer.Write(scanNum);
                writer.Write("\t");
                writer.Write(sequence);
                writer.Write("\t");
                writer.Write(mod);
                writer.Write("\t");
                writer.Write(mass);
                writer.Write("\t");
                writer.Write(protName);
                writer.Write("\t");
                writer.Write(protDesc);
                writer.Write("\t");
                writer.Write(firstResId);
                writer.Write("\t");
                writer.Write(lastResId);
                writer.Write("\t");
                writer.Write(score);
                writer.Write("\t");
                writer.Write(qvalue);
                writer.Write("\n");
            }
        }

        private void OutputAlignmentResult(LcMsFeatureAlignment align, string outFilePath, List<string> rawFiles, bool isTemp = true)
        {
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(outFilePath);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++)
            {
                var dataSetName = Path.GetFileNameWithoutExtension(rawFiles[i]);
                writer.Write("\t{0}", dataSetName);
            }

            for (var i = 0; i < align.CountDatasets; i++)
            {
                //var dataSetName = Path.GetFileNameWithoutExtension(align.RawFileList[i]);
                writer.Write("\t{0}_Score", i);
            }

            /*
            for (var i = 0; i < align.CountDatasets; i++)
            {
                //var dataSetName = Path.GetFileNameWithoutExtension(align.RawFileList[i]);
                writer.Write("\t{0}_Net", i);
            }*/

            writer.Write("\n");
            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = GetMinMaxNet(features);

                writer.Write(@"{0}\t{1:0.00000}\t{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Score : 0d);
                }
                /*
                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    if (feature != null) writer.Write("{0:0.00000}", feature.MinNet);
                    else writer.Write(0);
                }

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    if (feature != null) writer.Write("{0:0.00000}", feature.MaxNet);
                    else writer.Write(0);
                }*/

                writer.Write("\n");
            }
            writer.Close();

            if (isTemp) return;

            var outDirectory = Path.GetDirectoryName(Path.GetFullPath(outFilePath));
            for (var i = 0; i < align.CountDatasets; i++)
            {
                var dataSetName = Path.GetFileNameWithoutExtension(rawFiles[i]);
                //writer.Write("\t{0}", dataSetName);
                // now output results!!
                var ms1ftFilePath = String.Format(@"{0}\{1}.aligned.ms1ft", outDirectory, dataSetName);
                var writer2 = new StreamWriter(ms1ftFilePath);
                writer2.WriteLine(LcMsFeatureFinderLauncher.GetHeaderString());

                for (var j = 0; j < align.CountAlignedFeatures; j++)
                {
                    var f1 = alignedFeatureList[j][i];
                    writer2.Write("{0}\t", j + 1);
                    writer2.WriteLine(LcMsFeatureFinderLauncher.GetString(f1));
                }
                writer2.Close();
            }
        }

        private void RunFeatureAlignment(List<string> ms1FtFiles, List<string> rawFiles, string outFilePath)
        {
            var runList = new List<LcMsRun>();

            foreach(var rawFile in rawFiles)
                runList.Add(new PbfLcMsRun(rawFile));

            var align = new LcMsFeatureAlignment(ms1FtFiles, runList, new LcMsFeatureAlignComparer(new Tolerance(10)));
            align.AlignFeatures();
            Console.WriteLine("# of aligned features = {0}", align.CountAlignedFeatures);
            var tempOutPath = outFilePath + ".tmp";
            OutputAlignmentResult(align, tempOutPath, rawFiles, true);

            align.RefineAbundance();
            OutputAlignmentResult(align, outFilePath, rawFiles, false);
        }

        private Tuple<int, int> GetMinMaxMs1ScanNum(LcMsRun run, double minTime, double maxTime)
        {
            var ms1ScanNums = run.GetMs1ScanVector();
            var minScanNum = -1;
            var maxScanNum = -1;

            for(var i = 1; i < ms1ScanNums.Length; i++)
            {
                var time = run.GetElutionTime(ms1ScanNums[i]);
                if (minScanNum < 0 && time > minTime)
                {
                    minScanNum = ms1ScanNums[i - 1];
                }
                if (maxScanNum < 0 && time > maxTime)
                {
                    maxScanNum = ms1ScanNums[i];
                    break;
                }
            }
            return new Tuple<int, int>(minScanNum, maxScanNum);
        }

        public static Tuple<double, int, double, double> GetMinMaxNet(IList<LcMsFeature> features)
        {
            //var minNet = 1.0d;
            //var maxNet = 0d;
            //var minElutionTime = double.MaxValue;
            //var maxElutionTime = 0d;
            var massList = new List<double>();
            var minElutionList = new List<double>();
            var maxElutionList = new List<double>();

            var charge = 0;
            foreach (var f in features)
            {
                if (f == null) continue;
                //minNet = Math.Min(minNet, f.MinNet);
                //maxNet = Math.Max(maxNet, f.MaxNet);
                //minElutionTime = Math.Min(minElutionTime, f.MinElutionTime);
                //maxElutionTime = Math.Max(maxElutionTime, f.MaxElutionTime);
                minElutionList.Add(f.MinElutionTime);
                maxElutionList.Add(f.MaxElutionTime);
                massList.Add(f.Mass);
                charge = f.RepresentativeCharge;
            }
            massList.Sort();
            var mass = massList[(int)(massList.Count * 0.5)];
            var minElutionTime = minElutionList[(int)(massList.Count * 0.5)];
            var maxElutionTime = maxElutionList[(int)(massList.Count * 0.5)];

            //return new Tuple<double, int, double, double>(mass, charge, minNet, maxNet);
            return new Tuple<double, int, double, double>(mass, charge, minElutionTime, maxElutionTime);
        }
    }
}
