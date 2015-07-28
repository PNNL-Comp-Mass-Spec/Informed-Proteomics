using System;
using System.Collections.Generic;
using System.Data.Odbc;
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
using NUnit.Framework;

namespace InformedProteomics.Test
{
    class TestLcMsFeatureAlignment
    {
        [Test]
        public void TestPeptidomics()
        {
            const string ms1ftFolder = @"D:\MassSpecFiles\Peptidome";
            const string rawFolder = @"\\protoapps\UserData\Jungkap\peptidome";
            string outFilePath = string.Format(@"{0}\aligned_features.tsv", ms1ftFolder);

            var fileEntries = Directory.GetFiles(ms1ftFolder);
            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();
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

        [Test]
        public void TestUtexLcMsFeatureAlign()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            const string ms1ftFolder = @"D:\MassSpecFiles\UTEX\ms1ft";

            const string quantOutFile = @"D:\MassSpecFiles\UTEX\features_promex.tsv";
            const string mspOutFile = @"D:\MassSpecFiles\UTEX\features_msp.tsv";
            const string msaOutFile = @"D:\MassSpecFiles\UTEX\features_msa.tsv";
            
            if (!Directory.Exists(rawFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, rawFolder);
                return;
            }

            var fileEntries = Directory.GetFiles(ms1ftFolder);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();

            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();
            //var prsmContainer = new PrSmContainer[dataset.Count];

            var mspMap = new ProteinSpectrumMathMap[dataset.Count];
            var msaMap = new ProteinSpectrumMathMap[dataset.Count];

            for (var i = 0; i < dataset.Count; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, dataset[i]);

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

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                rawFiles.Add(rawFile);
                ms1FtFiles.Add(ms1File);
                //prsmContainer[i] = new PrSmContainer(i, dataset[i]);
                mspMap[i] = new ProteinSpectrumMathMap(run, i, dataset[i]);
                msaMap[i] = new ProteinSpectrumMathMap(run, i, dataset[i]);

                // load identification results
                Console.WriteLine(dataset[i]);

                var path = string.Format(@"D:\MassSpecFiles\UTEX\MSPF\{0}_IcTda.tsv", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                mspMap[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                Console.WriteLine("\t[MSP] Total prsm = {0}, total unique proteins = {1}", mspMap[i].CountIdentifiedScans(), mspMap[i].CountIdentifiedUniqueProteoforms());

                path = string.Format(@"D:\MassSpecFiles\UTEX\MSA\{0}_MSAlign_ResultTable.txt", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                msaMap[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);
                Console.WriteLine("\t[MSA] Total prsm = {0}, total unique proteins = {1}", msaMap[i].CountIdentifiedScans(), msaMap[i].CountIdentifiedUniqueProteoforms());
                
            }

            if (rawFiles.Count == 0)
            {
                Console.WriteLine(@"Warning: No files were found in method {0}", methodName);
                return;
            }

            var align = new LcMsFeatureAlignment(ms1FtFiles, rawFiles);
            align.AlignFeatures();
            Console.WriteLine("{0} alignments ", align.CountAlignedFeatures);

            //align.RefineAbundance();
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(quantOutFile);

            var mspWriter = new StreamWriter(mspOutFile);
            var msaWriter = new StreamWriter(msaOutFile);
            
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}", dataset[i]);

            mspWriter.Write("FeatureID");
            msaWriter.Write("FeatureID");
            for (var i = 0; i < align.CountDatasets; i++)
            {
                mspWriter.Write("\tProteinName_{0}", i);
                mspWriter.Write("\tSequence_{0}", i);

                msaWriter.Write("\tProteinName_{0}", i);
                msaWriter.Write("\tSequence_{0}", i);
            }
            writer.Write("\n");
            msaWriter.Write("\n");
            mspWriter.Write("\n");

            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = GetMinMaxNet(features);
                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);
                
                mspWriter.Write(i + 1);
                msaWriter.Write(i + 1);

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);

                    ProteinSpectrumMatch prsm = null;
                    if (feature != null)
                    {
                        prsm = mspMap[j].FindByFeature(feature, new Tolerance(10));
                    }
                    else
                    {
                        prsm = mspMap[j].FindByFeature(minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4, new Tolerance(10));
                    }
                    mspWriter.Write("\t");
                    mspWriter.Write((prsm != null) ? prsm.ProteinName : "");
                    mspWriter.Write("\t");
                    mspWriter.Write((prsm != null) ? prsm.SequenceText : "");

                    prsm = null;
                    if (feature != null)
                    {
                        prsm = msaMap[j].FindByFeature(feature, new Tolerance(10));
                    }
                    else
                    {
                        prsm = msaMap[j].FindByFeature(minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4, new Tolerance(10));
                    }
                    msaWriter.Write("\t");
                    msaWriter.Write((prsm  != null) ? prsm.ProteinName : "");
                    msaWriter.Write("\t");
                    msaWriter.Write((prsm  != null) ? prsm.SequenceText : "");
                }

                mspWriter.Write("\n");
                msaWriter.Write("\n");
                writer.Write("\n");
            }
            writer.Close();
            mspWriter.Close();
            msaWriter.Close();
        }

        [Test]
        public void TestSpikeInLcMsFeatureAlign()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFolder = @"D:\MassSpecFiles\CPTAC_spike_in\raw";
            const string ms1ftFolder = @"D:\MassSpecFiles\CPTAC_spike_in\ms1ft";
            const string outFile = @"D:\MassSpecFiles\CPTAC_spike_in\aligned_features.tsv";

            /*const string fastaFilePath = @"D:\MassSpecFiles\UTEX\db\ID_004352_779A168F.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }*/
            if (!Directory.Exists(rawFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, rawFolder);
                return;
            }

            var fileEntries = Directory.GetFiles(ms1ftFolder);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();

            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();
            var prsmContainer = new PrSmContainer[dataset.Count];
            for (var i = 0; i < dataset.Count; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, dataset[i]);

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
                prsmContainer[i] = new PrSmContainer(i, dataset[i], 20);

                // load identification results
                Console.WriteLine(dataset[i]);

                var path = string.Format(@"D:\MassSpecFiles\CPTAC_spike_in\MSPF\{0}_IcTda.tsv", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                /*
                path = string.Format(@"D:\MassSpecFiles\UTEX\MSA\{0}_MSAlign_ResultTable.txt", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);
                */
                Console.WriteLine("Total prsm = {0}, total unique proteins = {1}", prsmContainer[i].CountIdentifiedScans(), prsmContainer[i].CountIdentifiedUniqueProteoforms());
            }

            if (rawFiles.Count == 0)
            {
                Console.WriteLine(@"Warning: No files were found in method {0}", methodName);
                return;
            }

            var align = new LcMsFeatureAlignment(ms1FtFiles, rawFiles);
            align.AlignFeatures();
            Console.WriteLine("{0} alignments ", align.CountAlignedFeatures);

            align.RefineAbundance();
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(outFile);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}", dataset[i]);

            writer.Write("\tMSPathFinder_Protein");
            writer.Write("\tMSPathFinder_Sequence");
            writer.Write("\tMSPathFinder_Score");

            writer.Write("\tMSAlign_Protein");
            writer.Write("\tMSAlign_Sequence");
            writer.Write("\tMSAlign_Score");

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
                    writer.Write("\t");
                    writer.Write((prsmRet[k] != null) ? prsmRet[k].Score : 0);
                }
                writer.Write("\n");
            }
            writer.Close();
        }


        [Test]
        public void TestCompRefLcMsFeatureAlign()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string dataFolder = @"D:\MassSpecFiles\CompRef";
            const string outFile = @"D:\MassSpecFiles\CompRef\aligned_features_new.tsv";
            const string fastaFilePath = @"D:\MassSpecFiles\CompRef\db\ID_003278_4B4B3CB1.fasta";


            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }
            if (!Directory.Exists(dataFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, dataFolder);
                return;
            }

            var fileEntries = Directory.GetFiles(dataFolder);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("pbf") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();

            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();
            var prsmContainer = new PrSmContainer[dataset.Count];
            for (var i = 0; i < dataset.Count; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", dataFolder, dataset[i]);
                var ms1File = string.Format(@"{0}\ms1ft\{1}.ms1ft", dataFolder, dataset[i]);

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
                
                prsmContainer[i] = new PrSmContainer(i, dataset[i]);

                // load identification results
                Console.WriteLine(dataset[i]);

                var path = string.Format(@"D:\MassSpecFiles\CompRef\MsPathFinderMerged\{0}_IcTda.tsv", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsPathFinder);
               
                path = string.Format(@"D:\MassSpecFiles\CompRef\MsAlign\{0}_MSAlign_ResultTable.txt", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);

                Console.WriteLine("Total prsm = {0}, total unique proteins = {1}", prsmContainer[i].CountIdentifiedScans(), prsmContainer[i].CountIdentifiedUniqueProteoforms());
            }

            if (rawFiles.Count == 0)
            {
                Console.WriteLine(@"Warning: No files were found in method {0}", methodName);
                return;
            }
            
            var align = new LcMsFeatureAlignment(ms1FtFiles, rawFiles);
            align.AlignFeatures();
            Console.WriteLine("{0} alignments ", align.CountAlignedFeatures);

            align.RefineAbundance();
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(outFile);
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
            writer.Close();
        }

        private void OutputAlignmentResult(LcMsFeatureAlignment align, string outFilePath)
        {
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(outFilePath);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++)
            {
                var dataSetName = Path.GetFileNameWithoutExtension(align.FeatureFileList[i]);
                writer.Write("\t{0}", dataSetName);
            }
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

                writer.Write("\n");
            }
            writer.Close();               
        }

        private void RunFeatureAlignment(List<string> ms1FtFiles, List<string> rawFiles, string outFilePath)
        {
            var align = new LcMsFeatureAlignment(ms1FtFiles, rawFiles);
            align.AlignFeatures();
            Console.WriteLine("# of aligned features = {0}", align.CountAlignedFeatures);
            var tempOutPath = outFilePath + ".bak";
            OutputAlignmentResult(align, tempOutPath);
            
            align.RefineAbundance();
            OutputAlignmentResult(align, outFilePath);
            //align.TryFillMissingFeature(alignedFeatureList);
        }


        private ProteinSpectrumMatch[] FindProteinSpectrumMatch(PrSmContainer prsmContainer, LcMsFeature[] alignedFeatures, LcMsRun run, double mass, double minTime, double maxTime)
        {
            var i = prsmContainer.DataId;
            var ret = new ProteinSpectrumMatch[2];
            List<ProteinSpectrumMatch> allPrsmList = null;

            if (alignedFeatures[i] != null)
            {
                allPrsmList = prsmContainer.FindByFeature(alignedFeatures[i], new Tolerance(10));
                if (allPrsmList.Count < 1)
                {
                    allPrsmList = prsmContainer.FindByFeature(alignedFeatures[i], new Tolerance(13));
                }

            }
            else
            {
                var scanRange = GetMinMaxMs1ScanNum(run, minTime, maxTime);
                allPrsmList = prsmContainer.FindByFeature(mass, scanRange.Item1, scanRange.Item2, new Tolerance(10));
                if (allPrsmList.Count < 1)
                {
                    allPrsmList = prsmContainer.FindByFeature(mass, scanRange.Item1, scanRange.Item2, new Tolerance(13));
                }
            }

            if (allPrsmList == null) return ret;

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

        private ProteinSpectrumMatch[] FindProteinSpectrumMatch(PrSmContainer[] prsmContainer, LcMsFeature[] alignedFeatures)
        {
            var ret = new ProteinSpectrumMatch[2];
            var allPrsmList = new List<ProteinSpectrumMatch>();
            for(var i = 0; i < prsmContainer.Length; i++)
            {
                if (alignedFeatures[i] == null) continue;

                var prsmList = prsmContainer[i].FindByFeature(alignedFeatures[i], new Tolerance(10));
                if (prsmList.Count < 1)
                {
                    prsmList = prsmContainer[i].FindByFeature(alignedFeatures[i], new Tolerance(13));
                }

                foreach (var prsm in prsmList)
                {
                    if(!allPrsmList.Contains(prsm)) allPrsmList.Add(prsm);
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

        private Tuple<double, int, double, double> GetMinMaxNet(IList<LcMsFeature> features)
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

    }
}
