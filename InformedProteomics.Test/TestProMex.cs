using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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
    public class TestProMex
    {
        [Test]
        public void TestQuantBottomUpData()
        {
            var m = 500;
            var theoreticalEnvelope = new IsotopeList(m, 30, 0.1);
        }
        


        [Test]
        public void TestMs1EvidenceScore()
        {
            const string TestRawFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01.pbf";
            const string TestResultFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01_IcTda.tsv";

            var run = PbfLcMsRun.GetLcMsRun(TestRawFile);
            var tsvParser = new TsvFileParser(TestResultFile);
            var featureFinder = new Ms2FeatureQuntification(run);

            for (var i = 0; i < tsvParser.NumData; i++)
            {
                var scan = int.Parse(tsvParser.GetData("Scan")[i]);
                var charge = int.Parse(tsvParser.GetData("Charge")[i]);
                var mass = double.Parse(tsvParser.GetData("Mass")[i]);
                var qvalue = double.Parse(tsvParser.GetData("QValue")[i]);

                var ms2Feature = new Ms2Feature()
                {
                    Charge = charge,
                    ScanNum = scan,
                    Mass = mass,
                };
                
                var score = featureFinder.GetMs1EvidenceScore(ms2Feature);
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", scan, mass, charge, qvalue, score);
            }   
        }
            


        
        [Test]
        public void TestPredictPTMfromMs1ft()
        {
            const string resultFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            const string ms1ftFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";

            var parser = new TsvFileParser(resultFilePath);
            var sequences = parser.GetData("Sequence");
            var modifications = parser.GetData("Modifications");
            var scanNums = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var qValues = parser.GetData("QValue").Select(s => Convert.ToDouble(s)).ToArray();
            var nMacthed = parser.GetData("#MatchedFragments");
            var aaSet = new AminoAcidSet();
            var ptmList = new List<Tuple<int, double, double>>();

            for (var i = 0; i < parser.NumData; i++)
            {
                if (qValues[i] > 0.01) continue;
                //var sequenceComp = aaSet.GetComposition(sequences[i]) + Composition.H2O;
                var seq = new Sequence(sequences[i], aaSet);
                var sequenceComp = seq.Composition + Composition.H2O;

                var modComposition = Composition.Zero;
                var modsStr = modifications[i];
                if (modsStr.Length == 0) continue;
                var mods = modsStr.Split(',');
                foreach (var modStr in mods.Where(str => str.Length > 0))
                {
                    var modName = modStr.Split()[0];
                    var mod = Modification.Get(modName);
                    modComposition += mod.Composition;
                }
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", scanNums[i], sequenceComp.Mass, modComposition.Mass, nMacthed[i], sequences[i], modsStr);
                //var compFromSeqAndMods = sequenceComp + modComposition;
                //Assert.True(compFromSeqAndMods.Equals(compositions[i]));

                ptmList.Add(new Tuple<int, double, double>(scanNums[i], sequenceComp.Mass, modComposition.Mass));
            }

            //var featureParser = new TsvFileParser(ms1ftFilePath);
            //var minScan = featureParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            //var maxScan = featureParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            //var monoMass = featureParser.GetData("MonoMass").Select(s => Convert.ToDouble(s)).ToArray();
        }
        
        [Test]
        public void TestGeneratingMs1FeatureFile()
        {
            //const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";
            const string specFilePath = @"D:\MassSpecFiles\test\QC_Shew_Intact_4_01Jan15_Bane_C2-14-08-02RZ.raw";
            
            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 3100;
            const int maxThreads = 10;
            
            var param = new Ms1FeatureFinderInputParameter
            {
                InputPath = specFilePath,
                MinSearchMass = minScanMass,
                MaxSearchMass = maxScanMass,
                MinSearchCharge = minScanCharge,
                MaxSearchCharge = maxScanCharge,
                MaxThreads = maxThreads
            };
            var featureFinder = new Ms1FeatureFinderLauncher(param);
            featureFinder.Run();
        }

        [Test]
        public void TestProMexFilter()
        {
            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            const string ms1FtPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath, 0.15);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }
        /*
        [Test]
        public void TestAlignProMexResults()
        {
            //var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            //var featureDir = @"D:\Test\Quant\Histone";
            //var outFile = @"D:\Test\Quant\Histone\aligned_features.tsv";
            //var dataset = new string[2] {"MZ20150405FG_MT_OXI", "MZ20150405FG_WT_OXI"};
            
            //var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1";
            //var featureDir = @"D:\Test\Quant\cptac_10rep";
            //var outFile = @"D:\Test\CPTAC_Intact_rep10.tsv";

            //var featureDir = @"D:\Test\Quant\spike_in";
            //var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\";
            //var outFile = @"D:\Test\Quant\spike_in\aligned_features.tsv";
            //var dataset = new string[5];
            //dataset[0] = "CPTAC_Intact_1to2_8Apr15_Bane_14-09-03RZa";
            //dataset[1] = "CPTAC_Intact_1to5_8Apr15_Bane_14-09-03RZ";
            //dataset[2] = "CPTAC_Intact_1to10_8Apr15_Bane_14-09-03RZ";
            //dataset[3] = "CPTAC_Intact_1to20_8Apr15_Bane_14-09-03RZ";
            //dataset[4] = "CPTAC_Intact_SpikeTest_8Apr15_Bane_14-09-03RZ";
            //var minDatasetIndex = 0;
            //var maxDatasetIndex = 4;
            
            //var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2014_3";
            //var featureDir = @"D:\Test\Quant\lewy";
            //var outFile = @"D:\Test\Quant\lewy\Lewy_Intact_by_ProMex.tsv";
            
            var rawDir = @"\\protoapps\UserData\Jungkap\SBEP";
            var featureDir = @"D:\Test\Quant\SK";
            var outFile = @"D:\Test\Quant\SK\SBEP_features.tsv";
            var dataset = new string[8];
            dataset[0] = "SBEP_STM_001_02222012_Aragon";
            dataset[1] = "SBEP_STM_001_02272012_Aragon";
            dataset[2] = "SBEP_STM_002_02222012_Aragon";
            dataset[3] = "SBEP_STM_002_02222012_Aragon_120226215929";
            dataset[4] = "SBEP_STM_003_02222012_Aragon";
            dataset[5] = "SBEP_STM_003_02222012_Aragon_120227062113";
            dataset[6] = "SBEP_STM_004_02222012_Aragon";
            dataset[7] = "SBEP_STM_004_02272012_Aragon";

            var minDatasetIndex = 0;
            var maxDatasetIndex = dataset.Length - 1;
            

            var rawFiles = new List<string>();
            var ms1ftFiles = new List<string>();

            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++)
            {
                //var rawFile = string.Format(@"{0}\CPTAC_Intact_rep{1}_15Jan15_Bane_C2-14-08-02RZ.pbf", rawDir, i);
                //var ms1File = string.Format(@"{0}\CPTAC_Intact_rep{1}_15Jan15_Bane_C2-14-08-02RZ.ms1ft", featureDir, i);
                //var rawFile = string.Format(@"{0}\Lewy_intact_{1}.pbf", rawDir, i.ToString("D2"));
                //var ms1File = string.Format(@"{0}\Lewy_intact_{1}.ms1ft", featureDir, i.ToString("D2"));
                var rawFile = string.Format(@"{0}\{1}.pbf", rawDir, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", featureDir, dataset[i]);

                rawFiles.Add(rawFile);
                ms1ftFiles.Add(ms1File);
            }
            
            var align = new Ms1FeatureAlign(ms1ftFiles, rawFiles);
            align.AlignNetAll();

            var nDataset = maxDatasetIndex - minDatasetIndex + 1;
            var fid2gid = new int[nDataset][];
            for (var i = 0; i < nDataset; i++)
            {
                fid2gid[i] = new int[align.GetFeatureCount(i)];
            }

            double[,] abundanceData;
            var alignedFeatureList = align.ClusterFeatures(out abundanceData);
            
            var writer = new StreamWriter(outFile);

            writer.Write("MonoMass\tMinNet\tMaxNet");
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}", dataset[i]);
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}_featureId", dataset[i]);
            writer.Write("\n");

            for (var i = 0; i < alignedFeatureList.Count(); i++)
            {
                writer.Write(@"{0}	{1:0.00000}	{2:0.00000}", alignedFeatureList[i].RepresentativeMass, alignedFeatureList[i].MinNet, alignedFeatureList[i].MaxNet);
                //writer.Write("\t{0}", ArrayUtil.ToString(alignedFeatureList[i].GetMasses(), ";"));

                var sb = new StringBuilder();
                for (var j = 0; j < abundanceData.GetLength(1); j++)
                {
                    writer.Write("\t");
                    writer.Write(abundanceData[i, j]);
                }
                
                for (var j = 0; j < abundanceData.GetLength(1); j++)
                {
                    var featureIds = alignedFeatureList[i].GetFeatureId(j);
                    var fidList = ArrayUtil.ToString(featureIds, ";");
                    writer.Write("\t");
                    writer.Write(fidList);
                    foreach(var k in featureIds) fid2gid[j][k - 1] = i + 1;
                }

                writer.Write("\n");
            }
            writer.Close();

            for (var i = 0; i < nDataset; i++)
            {
                var wt = new StreamWriter(string.Format(@"{0}\{1}.fid2gid", featureDir, dataset[i]));
                for(var j = 0; j < fid2gid[i].Length; j++) wt.WriteLine(fid2gid[i][j]);
                wt.Close();
            }

            foreach (var et in align.ElutionLengths)
            {
                Console.WriteLine("{0}\t{1}", et.Item1, et.Item2);
            }

        }
         */
    }
    


}
