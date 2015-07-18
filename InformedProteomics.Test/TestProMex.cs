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
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;


namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestProMex
    {
        /*
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
       */ 

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
            var featureFinder = new LcMsPeakMatrix(run);

            for (var i = 0; i < tsvParser.NumData; i++)
            {
                var scan = int.Parse(tsvParser.GetData("Scan")[i]);
                var charge = int.Parse(tsvParser.GetData("Charge")[i]);
                var mass = double.Parse(tsvParser.GetData("Mass")[i]);
                var qvalue = double.Parse(tsvParser.GetData("QValue")[i]);

                //var targetFeature = new TargetFeature(mass, charge, scan);
                var score = featureFinder.GetMs1EvidenceScore(mass, scan, charge);
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

            const string outFolderPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output";
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

            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }

    }
    


}
