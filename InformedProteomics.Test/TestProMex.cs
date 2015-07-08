using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
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
    public class TestProMex
    {

        [Test]
        public void CollectTrainingSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string idFileFolder = @"\\protoapps\UserData\Jungkap\TrainingSet";
            if (!Directory.Exists(idFileFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
                return;
            }

            const string outFileFolder = @"\\protoapps\UserData\Jungkap\TrainingSet";
            var rawFileLists = new string[]
            {
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf"
            };

            //var trainingSetList = File.ReadAllLines(@"D:\MassSpecFiles\training\training_datasets.txt");
            //var trainSet = trainingSetList.Where(set => !set.StartsWith("#")).ToList();

            foreach (var dataset in rawFileLists)
            {
                if (!File.Exists(dataset))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", dataset);
                    continue;
                }

                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var idFile = string.Format(@"{0}\{1}_IcTda.tsv", idFileFolder, dataname);

                Console.WriteLine(dataset);
                Console.WriteLine(idFile);
                var targetSets = TargetFeature.CollectUniqueTargets(dataset, idFile, 0.005);
                Console.WriteLine(targetSets.Count);

                var writer = new StreamWriter(string.Format(@"{0}\{1}.trainset.tsv", outFileFolder, Path.GetFileNameWithoutExtension(dataset)));
                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMass\tSequence\tModifications\tComposition");

                foreach (var feature in targetSets)
                {
                    if (feature.Rows.Count < 2) continue;

                    writer.Write(feature.MinScanNum); writer.Write("\t");
                    writer.Write(feature.MaxScanNum); writer.Write("\t");
                    writer.Write(feature.MinCharge); writer.Write("\t");
                    writer.Write(feature.MaxCharge); writer.Write("\t");
                    writer.Write(feature.Mass); writer.Write("\t");
                    writer.Write(feature.Sequence); writer.Write("\t");
                    writer.Write(feature.Modifications); writer.Write("\t");
                    writer.Write(feature.Composition); writer.Write("\n");
                }

                writer.Close();
            }
        }
        
        [Test]
        public void TestQuantTopDownData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFile = @"D:\Test\Quant\spike_in\raw\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ.pbf";
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFile);

            //var target = new TargetFeature(8445.573787, 9, 3084);
            var target = new TargetFeature(5588.966676, 7, 3442);
            
            var quantAnalyzer = new TargetMs1FeatureMatrix(run);

            var ms1Feature = quantAnalyzer.FindMs1Feature(target);

            Console.WriteLine("{0}-{1}", ms1Feature.MinCharge, ms1Feature.MaxCharge);
            Console.WriteLine("{0}-{1}", ms1Feature.MinScanNum, ms1Feature.MaxScanNum);
            Console.WriteLine("{0}", ms1Feature.Abundance);

        }
        
        [Test]
        public void TestQuantBottomUpData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            for (var m = 500; m < 5000; m += 100)
            {
                var theoreticalEnvelope = new IsotopeList(m, 30, 0.1);    

                Console.WriteLine("{0}\t{1}", m, theoreticalEnvelope.Count);
            }
            
        }

        [Test]
        public void TestMs1EvidenceScore()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string TestRawFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01.pbf";
            if (!File.Exists(TestRawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestRawFile);
                return;
            }

            const string TestResultFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01_IcTda.tsv";
            if (!File.Exists(TestResultFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestResultFile);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(TestRawFile);
            var tsvParser = new TsvFileParser(TestResultFile);
            var featureFinder = new TargetMs1FeatureMatrix(run);

            for (var i = 0; i < tsvParser.NumData; i++)
            {
                var scan = int.Parse(tsvParser.GetData("Scan")[i]);
                var charge = int.Parse(tsvParser.GetData("Charge")[i]);
                var mass = double.Parse(tsvParser.GetData("Mass")[i]);
                var qvalue = double.Parse(tsvParser.GetData("QValue")[i]);

                var targetFeature = new TargetFeature(mass, charge, scan);
                
                var score = featureFinder.GetMs1EvidenceScore(targetFeature);
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", scan, mass, charge, qvalue, score);
            }   
        }
        
        [Test]
        public void TestPredictPTMfromMs1ft()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string resultFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            if (!File.Exists(resultFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, resultFilePath);
                return;
            }

            // const string ms1ftFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";

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
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\QC_ShewIntact_1_19Jun15_Bane_14-09-01RZ.pbf";
            if (!File.Exists(specFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, specFilePath);
                return;
            }

            const string outFolderPath = @"\\protoapps\UserData\Sangtae\TestData\Output";
            if (!Directory.Exists(outFolderPath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, outFolderPath);
                return;
            }

            //const string specFilePath = @"D:\MassSpecFiles\test\QC_Shew_Intact_4_01Jan15_Bane_C2-14-08-02RZ.raw";
            
            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 5000;
            const int maxThreads = 10;
            
            var param = new Ms1FeatureFinderInputParameter
            {
                InputPath = specFilePath,
                OutputPath = outFolderPath,
                MinSearchMass = minScanMass,
                MaxSearchMass = maxScanMass,
                MinSearchCharge = minScanCharge,
                MaxSearchCharge = maxScanCharge,
                CsvOutput = true,
                ScoreReport = false,
                MaxThreads = maxThreads
            };
            var featureFinder = new Ms1FeatureFinderLauncher(param);
            featureFinder.Run();
        }

        [Test]
        public void TestProMexFilter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            if (!File.Exists(specFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, specFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            
            const string ms1FtPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            if (!File.Exists(ms1FtPath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, ms1FtPath);
                return;
            }

            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath, 0.15);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }


        [Test]
        public void TestAlignProMexResults()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var featureDir = @"\\protoapps\UserData\Sangtae\TestData\Output";
            if (!Directory.Exists(featureDir))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, featureDir);
                return;
            }

            var rawDir = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2";
            var outFile = @"\\protoapps\UserData\Sangtae\TestData\Output\aligned_features.tsv";

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

            var minDatasetIndex = 0;
            var maxDatasetIndex = dataset.Count - 1;

            var rawFiles = new List<string>();
            var ms1ftFiles = new List<string>();

            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", rawDir, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", featureDir, dataset[i]);

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
                ms1ftFiles.Add(ms1File);

                Console.WriteLine(dataset[i]);
            }

            if (rawFiles.Count == 0)
            {
                Console.WriteLine(@"Warning: No files were found in method {0}", methodName);
                return;
            }

            var align = new Ms1FeatureAlign(ms1ftFiles, rawFiles);
            var alignedFeatureList = align.GroupFeatures();

            Console.WriteLine("{0} alignments ",alignedFeatureList.Count);
            
            var writer = new StreamWriter(outFile);
            writer.Write("MonoMass\tMinNet\tMaxNet");
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}", i);
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}_minScan", i);
            for (var i = minDatasetIndex; i <= maxDatasetIndex; i++) writer.Write("\t{0}_maxScan", i);

            writer.Write("\n");

            foreach(var featureSet in alignedFeatureList)
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
        public void ProMexDebug()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            for (var mass = 500; mass < 100000; mass += 500)
            {
                var isotopes = new IsotopeList(mass, 0.01d);
                //var isotopes = new IsotopeList(mass, 1000, 0.1);

                Console.Write(mass);
                Console.Write("\t");
                Console.Write(isotopes.Count);
                Console.Write("\t");

                //for (var charge = 2; charge <= 60; charge++)
                //{
                var charge = 10;


                var minMz = Ion.GetIsotopeMz(mass, charge, isotopes.First().Index);
                var maxMz = Ion.GetIsotopeMz(mass, charge, isotopes.Last().Index);
                var ppm = ((maxMz - minMz) / minMz) * 1e6;

                Console.Write(isotopes.First().Index);
                Console.Write("\t");
                Console.Write(isotopes.Last().Index);
                Console.Write("\t");
                Console.Write(minMz);
                Console.Write("\t");
                Console.Write(maxMz);
                Console.Write("\t");
                Console.Write(ppm);
                Console.Write("\t");
                //}


                //foreach(var p in isotopes.EnvelopePdf) Console.Write("{0}\t", p);

                Console.Write("\n");
            }

        }

    }
    


}
