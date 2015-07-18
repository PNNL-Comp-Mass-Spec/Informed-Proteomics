using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestUtex
    {
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
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, featureDir);
                return;
            }

            if (!Directory.Exists(mspDir))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, mspDir);
                return;
            }

            var dataset = GetDataList(featureDir);

            var tsvParser = new TsvFileParser(outFile);
            var massList = new List<double>();
            for (var i = 0; i < tsvParser.NumData; i++)
                massList.Add(double.Parse(tsvParser.GetData("MonoMass")[i]));

            var featureIdMap = new Dictionary<int, string>();
            var tolerance = new Tolerance(12);
            var headers = new List<string>();

            //foreach (var data in dataset)
            for (var d = 0; d < dataset.Count; d++)
            {
                var data = dataset[d];
                var minScanColName = string.Format("{0}_minScan", d);
                var maxScanColName = string.Format("{0}_maxScan", d);
                
                var fname = string.Format(@"{0}\{1}_IcTda.tsv", mspDir, data);
                var idParser = new TsvFileParser(fname);
                var idRows = idParser.GetRows();
                if (headers.Count < 1) headers.AddRange(idParser.GetHeaders());

                for (var i = 0; i < idParser.NumData; i++)
                {
                    var scan = int.Parse(idParser.GetData("Scan")[i]);
                    var mass = double.Parse(idParser.GetData("Mass")[i]);
                    var qvalue = double.Parse(idParser.GetData("QValue")[i]);

                    if (qvalue > 0.01) break;

                    var massTol = tolerance.GetToleranceAsTh(mass);

                    var idx = massList.BinarySearch(mass);
                    if (idx < 0) idx = ~idx;

                    var found = false;
                    for (var j = idx; j >= 0; j--)
                    {
                        if (Math.Abs(mass - massList[j]) > massTol) break;

                        if (tsvParser.GetData(minScanColName)[j].Length < 1) continue;

                        if (int.Parse(tsvParser.GetData(minScanColName)[j]) < scan && scan < int.Parse(tsvParser.GetData(maxScanColName)[j]))
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
                        if (int.Parse(tsvParser.GetData(minScanColName)[j]) < scan && scan < int.Parse(tsvParser.GetData(maxScanColName)[j]))
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
            writer.Write(ArrayUtil.ToString(headers.ToArray(), "\t"));
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
                    writer.Write("\t"); writer.Write("{0}", tsvParser.GetData(string.Format("{0}", i))[key]); 
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
            var dmsDir = @"\\proto-4\VOrbiETD02\2015_2";

            if (!Directory.Exists(featureDir))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, featureDir);
                return;
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
                var dname = string.Format(@"{0}\{1}", dmsDir, data);
                var destFname = string.Format(@"{0}\MSP\{1}.zip", featureDir, data);

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
