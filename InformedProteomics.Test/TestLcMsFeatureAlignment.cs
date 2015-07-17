using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    class TestLcMsFeatureAlignment
    {
        [Test]
        public void TestCompRefLcMsFeatureAlign()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string dataFolder = @"D:\MassSpecFiles\CompRef";
            const string outFile = @"D:\MassSpecFiles\CompRef\aligned_features_new.tsv";
            const string fastaFilePath = @"D:\MassSpecFiles\CompRef\ID_003278_4B4B3CB1.fasta";


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

            var dataset = (from fileName in fileEntries where fileName.EndsWith("ms1ft") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();

            var rawFiles = new List<string>();
            var ms1FtFiles = new List<string>();
            var prsmContainer = new PrSmContainer[dataset.Count];
            for (var i = 0; i < dataset.Count; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", dataFolder, dataset[i]);
                var ms1File = string.Format(@"{0}\ms1ft_new\{1}.ms1ft", dataFolder, dataset[i]);

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
                
                prsmContainer[i] = new PrSmContainer(i, dataset[i], fastaFilePath);

                // load identification results
                Console.WriteLine(dataset[i]);

                var path = string.Format(@"D:\MassSpecFiles\CompRef\ProMexV2\{0}_IcTda.tsv", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                prsmContainer[i].AddIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                
                path = string.Format(@"D:\MassSpecFiles\CompRef\seqtag\{0}_tagmatch.tsv", dataset[i]);
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
            var alignedFeatureList = align.GroupFeatures();
            Console.WriteLine("{0} alignments ", align.CountAlignedFeatures);

            //QuantifyAgain(ref align);

            var writer = new StreamWriter(outFile);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}", dataset[i]);

            writer.Write("\tProteinName1");
            writer.Write("\tSequence1");

            writer.Write("\tProteinName2");
            writer.Write("\tSequence2");

            writer.Write("\tProteinName3");
            writer.Write("\tSequence3");

            //for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}_Abu2", i);
            //for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}_Score", i);
            //for (var i = 0; i < align.CountDatasets; i++) writer.Write("\t{0}_maxScan", i);

            writer.Write("\n");
            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = GetMinMaxNet(features);

                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);

                //var sb = new StringBuilder();
                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }

                var prsmList = FindProteinSpectrumMatch(prsmContainer, features);

                for (var k = 0; k < 3; k++)
                {
                    if (k < prsmList.Count)
                    {
                        writer.Write("\t");
                        writer.Write(prsmList[k].ProteinName);
                        writer.Write("\t");
                        writer.Write(prsmList[k].SequenceText);
                    }
                    else
                    {
                        writer.Write("\t");
                        writer.Write("");
                        writer.Write("\t");
                        writer.Write("");
                    }
                }

                if (prsmList.Count > 0) Console.WriteLine(i+1);
                /*
                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    if (feature != null) writer.Write(feature.AbundanceForBestCharges);
                    else writer.Write("0");
                }
                */
                writer.Write("\n");
            }
            writer.Close();
            
        }

        private List<ProteinSpectrumMatch> FindProteinSpectrumMatch(PrSmContainer[] prsmContainer, LcMsFeature[] alignedFeatures)
        {
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

            if (allPrsmList.Count < 1) return allPrsmList;

            allPrsmList.Sort();
            //return allPrsmList[0];
            return allPrsmList;

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

        private void QuantifyAgain(ref LcMsFeatureAlignment alignment)
        {
            var nDataset = alignment.CountDatasets;
            var featureGroup = alignment.GroupFeatures();
            var nFeatures = alignment.CountAlignedFeatures;

            for (var i = 0; i < nDataset; i++)
            {
                var run = PbfLcMsRun.GetLcMsRun(alignment.RawFileList[i]);
                var ms1ScanNums = run.GetMs1ScanVector();
                var featureFinder = new LcMsPeakMatrix(run);

                for (var j = 0; j < nFeatures; j++)
                {
                    var minNet = 0d;
                    var maxNet = 0d;
                    var mass = 0d;
                    var charge = 0;

                    if (featureGroup[j][i] == null)
                    {
                        var netPair = GetMinMaxNet(featureGroup[j]);
                        mass = netPair.Item1;
                        charge = netPair.Item2;
                        minNet = netPair.Item3;
                        maxNet = netPair.Item4;
                    }
                    else
                    {
                        mass = featureGroup[j][i].Mass;
                        charge = featureGroup[j][i].RepresentativeCharge;
                        minNet = featureGroup[j][i].MinNet;
                        maxNet = featureGroup[j][i].MaxNet;
                    }

                    var minScanNum = -1;
                    var maxScanNum = ms1ScanNums.Last();
                    for (var k = 0; k < ms1ScanNums.Length; k++)
                    {
                        var net = run.GetElutionTime(ms1ScanNums[k]) / run.GetElutionTime(run.MaxLcScan);
                        if (net > minNet && minScanNum < 0)
                        {
                            minScanNum = (k == 0) ? ms1ScanNums[k] : ms1ScanNums[k - 1];
                        }

                        if (net > maxNet)
                        {
                            maxScanNum = ms1ScanNums[k];
                            break;
                        }
                    }
                    if (minScanNum < 0) minScanNum = 0;

                    var feature = featureFinder.GetLcMsPeakCluster(mass, charge, minScanNum, maxScanNum);

                    if (feature != null && feature.Score > -10) featureGroup[j][i] = feature;
                }

                Console.WriteLine("{0} has been processed...", alignment.RawFileList[i]);
            }
        }
    
    }
}
