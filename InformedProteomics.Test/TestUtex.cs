using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Quantification;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestUtex
    {
        [Test]
        public void TestTagAlignedFeatures()
        {
            var mspDir = @"D:\MassSpecFiles\UTEX\MSP";
            var outFile = @"D:\MassSpecFiles\UTEX\aligned_features.tsv";
            var resultFile = @"D:\MassSpecFiles\UTEX\aligned_ids.tsv";
            var dataset = GetDataList();

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


        private List<string> GetDataList()
        {
            var featureDir = @"D:\MassSpecFiles\UTEX";
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
        public void TestFeatureAlign()
        {
            var featureDir = @"D:\MassSpecFiles\UTEX";
            var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            var outFile = @"D:\MassSpecFiles\UTEX\aligned_features.tsv";

            var fileEntries = Directory.GetFiles(featureDir);

            var dataset = GetDataList();
            
            var minDatasetIndex = 0;
            var maxDatasetIndex = dataset.Count - 1;

            var rawFiles = new List<string>();
            var ms1ftFiles = new List<string>();

            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawDir, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", featureDir, dataset[i]);

                rawFiles.Add(rawFile);
                ms1ftFiles.Add(ms1File);

                Console.WriteLine(dataset[i]);
            }

            var align = new Ms1FeatureAlign(ms1ftFiles, rawFiles);
            var alignedFeatureList = align.GroupFeatures();

            Console.WriteLine("{0} alignments ", alignedFeatureList.Count);

            var writer = new StreamWriter(outFile);
            writer.Write("MonoMass\tMinNet\tMaxNet");
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}", i);
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}_minScan", i);
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}_maxScan", i);

            writer.Write("\n");

            foreach (var featureSet in alignedFeatureList)
            {
                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", featureSet[0].RepresentativeMass, featureSet[0].MinNet, featureSet[0].MaxNet);

                var features = new Ms1Feature[dataset.Count];
                foreach (var f in featureSet) features[f.DataSetId] = f;


                var sb = new StringBuilder();
                for (var j = minDatasetIndex; j <= maxDatasetIndex; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }

                for (var j = minDatasetIndex; j <= maxDatasetIndex; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    if (feature != null) writer.Write(feature.MinScanNum);
                    else writer.Write("");
                }

                for (var j = minDatasetIndex; j <= maxDatasetIndex; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    if (feature != null) writer.Write(feature.MaxScanNum);
                    else writer.Write("");
                }

                writer.Write("\n");
            }
            writer.Close();
        }
        
        [Test]
        public void CopyUTEX()
        {
            var featureDir = @"D:\MassSpecFiles\UTEX";
            var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            var outFile = @"D:\MassSpecFiles\UTEX\aligned_features.tsv";
            var dmsDir = @"\\proto-4\VOrbiETD02\2015_2";

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
