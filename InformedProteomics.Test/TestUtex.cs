using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using NUnit.Framework;
using ProMex;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestUtex
    {
        private const string ProteinNamePrefix = "M744_";

        [Test]
        public void TestQuantifyIdedProteoforms()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            const string promexOutFolder = @"D:\MassSpecFiles\UTEX\MSAlign";
            const string msAlignResultFolder = @"D:\MassSpecFiles\UTEX\MSAlign";

            if (!Directory.Exists(rawFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, rawFolder);
            }

            var nDataset = 32;
            var dataset = new string[nDataset];
            for (var i = 0; i < nDataset; i++)
            {
                dataset[i] = String.Format("Syn_utex2973_Top_{0,2:D2}_TopDown_7May15_Bane_14-09-01RZ", i + 1);
                //var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
            }

            var prsmReader = new ProteinSpectrumMatchReader(0.01);

            var tolerance = new Tolerance(10);
            for (var i = 0; i < dataset.Length; i++)
            {
                var rawFile = String.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }
                var run = PbfLcMsRun.GetLcMsRun(rawFile);

                var path = String.Format(@"{0}\{1}_MSAlign_ResultTable.txt", msAlignResultFolder, dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }

                var prsmList = prsmReader.LoadIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);

                for (var j = 0; j < prsmList.Count; j++)
                {
                    var match = prsmList[j];
                    match.ProteinId = match.ProteinName.Substring(match.ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
                }
                
                // PrSM To Feature
                var prsmToFeatureIdMap = new int[prsmList.Count];
                for (var k = 0; k < prsmToFeatureIdMap.Length; k++) prsmToFeatureIdMap[k] = -1;

                // Feature To PrSM
                var featureToPrsm = new List<ProteinSpectrumMatchSet>();

                var featureFinder = new LcMsPeakMatrix(run, new LcMsFeatureLikelihood());
                var featureList = new List<LcMsPeakCluster>();
                var featureId = 0;
                for(var j = 0; j < prsmList.Count; j++)
                {
                    if (prsmToFeatureIdMap[j] >= 0) continue;
                    
                    var match = prsmList[j];
                    var minScanNum = match.ScanNum;
                    var maxScanNum = match.ScanNum;
                    var mass = match.Mass;
                    var charge = match.Charge;
                    var massTh = tolerance.GetToleranceAsTh(mass);
                    var id1 = match.ProteinId;

                    var feature = featureFinder.GetLcMsPeakCluster(mass, charge, minScanNum, maxScanNum);
                    var prsmSet = new ProteinSpectrumMatchSet(i) { match };
                    if (feature == null)
                    {
                        feature = featureFinder.CollectLcMsPeaksWithNoise(mass, charge, minScanNum, maxScanNum, charge, charge);
                        prsmToFeatureIdMap[j] = featureId;
                    }
                    else
                    {
                        prsmToFeatureIdMap[j] = featureId;
                        var etTol = Math.Max(run.GetElutionTime(run.MaxLcScan)*0.005, feature.ElutionLength*0.2);

                        for (var k = j + 1; k < prsmList.Count; k++)
                        {
                            var otherMatch = prsmList[k];
                            var id2 = otherMatch.ProteinId;
                            var et2 = run.GetElutionTime(otherMatch.ScanNum);

                            if (id1.Equals(id2) &&
                                feature.MinElutionTime - etTol < et2 && et2 < feature.MaxElutionTime - etTol &&
                                Math.Abs(otherMatch.Mass - mass) < massTh)
                            {
                                prsmToFeatureIdMap[k] = featureId;
                                prsmSet.Add(otherMatch);
                            }
                        }
                        
                    }
                    featureId++;

                    feature.Flag = 1;
                    featureList.Add(feature);
                    featureToPrsm.Add(prsmSet);
                }

                // overalp between features???
                for (var j = 0; j < featureList.Count; j++)
                {
                    var f1 = featureList[j];
                    if (f1.Flag < 1) continue;
                    var prsm1 = featureToPrsm[j];

                    for (var k = j+1; k < featureList.Count; k++)
                    {
                        var f2 = featureList[k];
                        if (f2.Flag < 1) continue;

                        var prsm2 = featureToPrsm[k];
                        if (Math.Abs(f1.Mass - f2.Mass) > tolerance.GetToleranceAsTh(f1.Mass)) continue;
                        if (!f1.CoElutedByNet(f2, 0.005)) continue;
                        if (!prsm1.ShareProteinId(prsm2)) continue;

                        // let us merge!!
                        if (f1.ScanLength > f2.ScanLength)
                        {
                            prsm1.AddRange(prsm2);
                            prsm2.Clear();
                            f2.Flag = 0;
                        }
                        else
                        {
                            prsm2.AddRange(prsm1);
                            prsm1.Clear();
                            f1.Flag = 0;
                        }
                    }
                }

                // now output results!!                
                var ms1ftFilePath = String.Format(@"{0}\{1}.ms1ft", promexOutFolder, dataset[i]);
                var writer = new StreamWriter(ms1ftFilePath);
                writer.WriteLine(LcMsFeatureFinderLauncher.GetHeaderString());

                for (var j = 0; j < featureList.Count; j++)
                {
                    var f1 = featureList[j];
                    if (f1.Flag < 1) continue;
                    var prsm1 = featureToPrsm[j];

                    var minScanNum = run.GetPrevScanNum(prsm1.MinScanNum, 1);
                    var maxScanNum = run.GetNextScanNum(prsm1.MaxScanNum, 1);
                    f1.ExpandScanRange(minScanNum, maxScanNum);
                    
                    writer.Write("{0}\t", j+1);
                    writer.WriteLine(LcMsFeatureFinderLauncher.GetString(f1));
                }
                writer.Close();

                Console.WriteLine(ms1ftFilePath);
            }
        }

        [Test]
        public void TestAlignFeatures()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            const string promexOutFolder = @"D:\MassSpecFiles\UTEX\MSAlign";
            const string msAlignResultFolder = @"D:\MassSpecFiles\UTEX\MSAlign";

            if (!Directory.Exists(rawFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, rawFolder);
            }

            var nDataset =32;
            var dataset = new string[nDataset];
            for (var i = 0; i < nDataset; i++)
            {
                dataset[i] = String.Format("Syn_utex2973_Top_{0,2:D2}_TopDown_7May15_Bane_14-09-01RZ", i + 1);
                //var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
            }
            
            var tolerance = new Tolerance(10);
            var ftComparer = new UtexFeatureComparer(tolerance);
            var align = new LcMsFeatureAlignment(ftComparer);
            var prsmReader = new ProteinSpectrumMatchReader(0.01);
            var validCount = 0;

            for (var i = 0; i < dataset.Length; i++)
            {
                var rawFile = String.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }
                var run = PbfLcMsRun.GetLcMsRun(rawFile);

                var path = String.Format(@"{0}\{1}_MSAlign_ResultTable.txt", msAlignResultFolder, dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }

                var ms1ftPath = String.Format(@"{0}\{1}.ms1ft", promexOutFolder, dataset[i]);
                if (!File.Exists(ms1ftPath))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", ms1ftPath);
                    continue;
                }

                validCount++;

                //var map = new ProteinSpectrumMathMap(run, i, dataset[i]);
                //map.LoadIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);
                var prsmList = prsmReader.LoadIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);

                for (var j = 0; j < prsmList.Count; j++)
                {
                    var match = prsmList[j];
                    match.ProteinId =
                        match.ProteinName.Substring(
                            match.ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
                }


                var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1ftPath, run);
                
                // tag features by PrSMs
                for(var j = 0; j < features.Count; j++)
                {
                    //features[j].ProteinSpectrumMatches = new ProteinSpectrumMatchSet(i);
                    var massTol = tolerance.GetToleranceAsTh(features[j].Mass);
                    foreach (var match in prsmList)
                    {
                        if (features[j].MinScanNum < match.ScanNum && match.ScanNum < features[j].MaxScanNum && Math.Abs(features[j].Mass - match.Mass) < massTol)
                        {
                            features[j].ProteinSpectrumMatches.Add(match);
                        }
                    }
                }

                align.AddDataSet(i, features, run);
            }

            if (validCount == 0)
            {
                Assert.Ignore("No files found!");
            }

            align.AlignFeatures();
            Console.WriteLine("{0} alignments ", align.CountAlignedFeatures);
            align.RefineAbundance();

            var alignedFeatureList = align.GetAlignedFeatures();
            for (var i = 0; i < nDataset; i++)
            {
                var ms1ftPath = String.Format(@"{0}\{1}_aligned.ms1ft", promexOutFolder, dataset[i]);
                var writer = new StreamWriter(ms1ftPath);
                writer.Write(LcMsFeatureFinderLauncher.GetHeaderString());
                writer.WriteLine("\tIdedMs2ScanNums");

                for (var j = 0; j < alignedFeatureList.Count; j++)
                {
                    writer.Write(j + 1);
                    writer.Write("\t");

                    if (alignedFeatureList[j][i] == null)
                    {
                        for (var k = 0; k < 14; k++) writer.Write("0\t");
                        writer.Write("0\n");
                    }
                    else
                    {
                        writer.Write(LcMsFeatureFinderLauncher.GetString(alignedFeatureList[j][i]));
                        writer.Write("\t");

                        if (alignedFeatureList[j][i].ProteinSpectrumMatches == null)
                        {
                            writer.Write("");
                        }
                        else
                        {
                            var scanNums = string.Join(";", alignedFeatureList[j][i].ProteinSpectrumMatches.Select(prsm => prsm.ScanNum));
                            writer.Write(scanNums);
                        }
                        
                        writer.Write("\n");
                    }
                }
                writer.Close();
            }
            
        }

        internal class UtexFeatureComparer : INodeComparer<LcMsFeature>
        {
            public UtexFeatureComparer(Tolerance tolerance = null)
            {
                _tolerance = tolerance ?? new Tolerance(10);
            }

            public bool SameCluster(LcMsFeature f1, LcMsFeature f2)
            {
                if (f1.DataSetId == f2.DataSetId) return false;
                // tolerant in mass dimension?
                
                var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));
                if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;
                
                // tolerant in elution time dimension?
                //var lenDiff = Math.Abs(f1.NetLength - f2.NetLength) / Math.Min(f1.NetLength, f2.NetLength);
                //if (lenDiff > 0.5) return false;

                if (!f1.CoElutedByNet(f2, 0.004)) return false; //e.g) 200*0.001 = 0.2 min = 30 sec
                
                if (f1.ProteinSpectrumMatches.ShareProteinId(f2.ProteinSpectrumMatches)) return true;
                
                return false;
            }

            private readonly Tolerance _tolerance;
        }


        /*
       internal class MsAlignPrsmComparer : INodeComparer<ProteinSpectrumMatch>
       {
           internal MsAlignPrsmComparer(LcMsRun run)
           {
               _run = run;
               _elutionWindow = _run.GetElutionTime(_run.MaxLcScan)*0.005;
           }

           public bool SameCluster(ProteinSpectrumMatch node1, ProteinSpectrumMatch node2)
           {
               var id1 = node1.ProteinName.Substring(node1.ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
               var id2 = node2.ProteinName.Substring(node2.ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
               if (!id1.Equals(id2)) return false;

               var tolerance = new Tolerance(10);
               var massTh = tolerance.GetToleranceAsTh(node1.Mass);
               var massDiff = Math.Abs(node1.Mass - node2.Mass);
                
               if (massDiff > massTh) return false;

               var et1 = _run.GetElutionTime(node1.ScanNum);
               var et2 = _run.GetElutionTime(node2.ScanNum);

               if (Math.Abs(et1 - et2) > _elutionWindow) return false;

               return true;
           }

           private readonly double _elutionWindow;
           private readonly LcMsRun _run;
       }

       internal class MsAlignPrsmGroupComparer : INodeComparer<ProteinSpectrumMatcheSet>
       {
           internal MsAlignPrsmGroupComparer(LcMsRun[] run)
           {
               _runArray = run;
               _elutionWindow = _runArray[0].GetElutionTime(_runArray[0].MaxLcScan) * 0.005;
           }
            
           public bool SameCluster(ProteinSpectrumMatcheSet node1, ProteinSpectrumMatcheSet node2)
           {
               if (node1.DataId == node2.DataId) return false;
                
               var id1 = node1[0].ProteinName.Substring(node1[0].ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
               var id2 = node2[0].ProteinName.Substring(node2[0].ProteinName.IndexOf(ProteinNamePrefix) + ProteinNamePrefix.Length, 5);
               var tolerance = new Tolerance(10);
                
               if (!id1.Equals(id2)) return false;

               var massTh = tolerance.GetToleranceAsTh(node1[0].Mass);
               var massMatch = false;
               foreach (var n1 in node1)
               {
                   foreach (var n2 in node2)
                   {
                       var massDiff = Math.Abs(n1.Mass - n2.Mass);
                       //if (massDiff < massTh || Math.Abs(massDiff - 1) < massTh)
                       if (massDiff < massTh)
                       {
                           massMatch = true;
                           break;
                       }
                   }
                   if (massMatch) break;
               }
                
               if (!massMatch) return false;

               var et1Min = _runArray[node1.DataId].GetElutionTime(node1.MinScanNum);
               var et1Max = _runArray[node1.DataId].GetElutionTime(node1.MaxScanNum);
                
               var et2Min = _runArray[node2.DataId].GetElutionTime(node2.MinScanNum);
               var et2Max = _runArray[node2.DataId].GetElutionTime(node2.MaxScanNum);

               if (et1Min >= et2Min && et1Min <= et2Max) return true;
               if (et1Max >= et2Min && et1Max <= et2Max) return true;

               if (Math.Abs((et1Min + et1Max)*0.5 - (et2Min + et2Max)*0.5) < _elutionWindow) return true;
               if (Math.Abs(et1Min - et2Min) < _elutionWindow) return true;
               if (Math.Abs(et1Max - et2Max) < _elutionWindow) return true;
               if (Math.Abs(et1Min - et2Max) < _elutionWindow) return true;
               if (Math.Abs(et1Max - et2Min) < _elutionWindow) return true;

               return false;                
           }
           private const string ProteinNamePrefix = "M744_";
           private readonly double _elutionWindow;
           private readonly LcMsRun[] _runArray;
       }
        [Test]
        public void TestQunatifyMsAlignResult()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            const string msaOutFile = @"D:\MassSpecFiles\UTEX\msalign_ided.tsv";

            if (!Directory.Exists(rawFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, rawFolder);
            }

            var nDataset = 32;

            var dataset = new string[nDataset];
            var runArray = new LcMsRun[nDataset];
            for (var i = 0; i < nDataset; i++)
            {
                dataset[i] = String.Format("Syn_utex2973_Top_{0,2:D2}_TopDown_7May15_Bane_14-09-01RZ", i+1);
                var rawFile = String.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                runArray[i] = PbfLcMsRun.GetLcMsRun(rawFile);
            }

            //var msaMap = new ProteinSpectrumMathMap[dataset.Length];
            var alignement = new ProteinSpectrumMatchAlignment();
            var prsmGroup = new List<ProteinSpectrumMatcheSet>[runArray.Length];
            for (var i = 0; i < dataset.Length; i++)
            {
                var rawFile = String.Format(@"{0}\{1}.pbf", rawFolder, dataset[i]);
                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }

                var run = runArray[i];
                var map = new ProteinSpectrumMathMap(run, i, dataset[i]);
                // load identification results
                Console.WriteLine(dataset[i]);
                var path = String.Format(@"D:\MassSpecFiles\UTEX\MSA\{0}_MSAlign_ResultTable.txt", dataset[i]);
                if (!File.Exists(path))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", path);
                    continue;
                }
                
                map.LoadIdentificationResult(path, ProteinSpectrumMatch.SearchTool.MsAlign);
                

                var comparer = new MsAlignPrsmComparer(runArray[i]);
                prsmGroup[i] = alignement.GroupingByPrsm(i, map.ProteinSpectrumMatches, comparer);


                Console.WriteLine("\t[MSA] Total prsm = {0}, total unique proteins = {1}, grouped PRSMs = {2}", map.CountIdentifiedScans(), map.CountIdentifiedUniqueProteoforms(),
                    prsmGroup[i].Count);
            }
            
            var alignGroupComparer = new MsAlignPrsmGroupComparer(runArray);
            var alignedPrsmSet = alignement.GroupAcrossRuns(prsmGroup, alignGroupComparer);


            Console.WriteLine(alignedPrsmSet.Length);
        }
        */

        [Test]
        public void TestTagAlignedFeatures()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var featureDir = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output";
            var mspDir = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output\MSP";
            var outFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output\aligned_features.tsv";
            var resultFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output\aligned_ids.tsv";

            if (!Directory.Exists(featureDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, featureDir);
            }

            if (!Directory.Exists(mspDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, mspDir);
            }

            var dataset = GetDataList(featureDir);

            var tsvParser = new TsvFileParser(outFile);
            var massList = new List<double>();
            for (var i = 0; i < tsvParser.NumData; i++)
                massList.Add(Double.Parse(tsvParser.GetData("MonoMass")[i]));

            var featureIdMap = new Dictionary<int, string>();
            var tolerance = new Tolerance(12);
            var headers = new List<string>();

            //foreach (var data in dataset)
            for (var d = 0; d < dataset.Count; d++)
            {
                var data = dataset[d];
                var minScanColName = String.Format("{0}_minScan", d);
                var maxScanColName = String.Format("{0}_maxScan", d);
                
                var fname = String.Format(@"{0}\{1}_IcTda.tsv", mspDir, data);
                var idParser = new TsvFileParser(fname);
                var idRows = idParser.GetRows();
                if (headers.Count < 1) headers.AddRange(idParser.GetHeaders());

                for (var i = 0; i < idParser.NumData; i++)
                {
                    var scan = Int32.Parse(idParser.GetData("Scan")[i]);
                    var mass = Double.Parse(idParser.GetData("Mass")[i]);
                    var qvalue = Double.Parse(idParser.GetData("QValue")[i]);

                    if (qvalue > 0.01) break;

                    var massTol = tolerance.GetToleranceAsTh(mass);

                    var idx = massList.BinarySearch(mass);
                    if (idx < 0) idx = ~idx;

                    var found = false;
                    for (var j = idx; j >= 0; j--)
                    {
                        if (Math.Abs(mass - massList[j]) > massTol) break;

                        if (tsvParser.GetData(minScanColName)[j].Length < 1) continue;

                        if (Int32.Parse(tsvParser.GetData(minScanColName)[j]) < scan && scan < Int32.Parse(tsvParser.GetData(maxScanColName)[j]))
                        {
                            found = true;
                            if (!featureIdMap.ContainsKey(j)) featureIdMap.Add(j, idRows[i]);
                            break;
                        }
                    }
                    
                    if (found) continue;
                    for (var j = idx + 1; j < massList.Count; j++)
                    {
                        if (Math.Abs(mass - massList[j]) > massTol) break;
                        if (tsvParser.GetData(minScanColName)[j].Length < 1) continue;
                        if (Int32.Parse(tsvParser.GetData(minScanColName)[j]) < scan && scan < Int32.Parse(tsvParser.GetData(maxScanColName)[j]))
                        {
                            found = true;
                            if (!featureIdMap.ContainsKey(j)) featureIdMap.Add(j, idRows[i]);
                            break;
                        }
                    }
                }
            }
            
            var writer = new StreamWriter(resultFile);

            writer.Write("AlignedFeatureID"); writer.Write("\t");
            writer.Write(string.Join("\t", headers));
            for (var i = 0; i < 32; i++)
            {
                writer.Write("\t");  writer.Write("{0}", i); 
            }
            writer.Write("\n");

            var id = 1;
            foreach (var key in featureIdMap.Keys)
            {
                writer.Write(id); writer.Write("\t");
                writer.Write(featureIdMap[key]);
                for (var i = 0; i < 32; i++)
                {
                    writer.Write("\t"); writer.Write("{0}", tsvParser.GetData(String.Format("{0}", i))[key]); 
                }
                writer.Write("\n");
                id++;
            }
            writer.Close();

        }


        private List<string> GetDataList(string featureDir)
        {            
            var fileEntries = Directory.GetFiles(featureDir);
            var dataset = new List<string>();
            foreach (string fileName in fileEntries)
            {
                if (fileName.EndsWith("ms1ft"))
                {
                    dataset.Add(Path.GetFileNameWithoutExtension(fileName));
                }
            }
            dataset.Sort();
            return dataset;
        }
            
        [Test]
        public void CopyUTEX()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var featureDir = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output";
            //var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            //var outFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output\aligned_features.tsv";
            var dmsDir = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\ProMex";

            if (!Directory.Exists(featureDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, featureDir);
            }

            var fileEntries = Directory.GetFiles(featureDir);

            var dataset = new List<string>();
            foreach (string fileName in fileEntries)
            {
                if (fileName.EndsWith("ms1ft"))
                {
                    dataset.Add(Path.GetFileNameWithoutExtension(fileName));

                }
            }
            dataset.Sort();

            foreach (var data in dataset)
            {
                var dname = String.Format(@"{0}\{1}", dmsDir, data);
                var destFname = String.Format(@"{0}\MSP\{1}.zip", featureDir, data);

                if (File.Exists(destFname))
                {
                    File.Delete(destFname);
                }

                var directories = Directory.GetDirectories(dname);
                //Console.WriteLine(dname);
                foreach (var dir in directories)
                {
                    var subDir = Path.GetFileName(dir);

                    if (subDir.StartsWith("MSP"))
                    {
                        //var dname2 = string.Format(@"{0}\{1}\{2}", dmsDir, data, dir);
                        var mspFiles = Directory.GetFiles(dir);

                        foreach (var fname in mspFiles)
                        {
                            if (fname.EndsWith("zip"))
                            {
                                Console.WriteLine(fname);
                                File.Copy(fname, destFname);
                            }
                        }

                    }
                }

            }
        }
    }
}
