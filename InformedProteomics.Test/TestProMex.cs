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
using InformedProteomics.Graphics;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;
using ProMex;


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

            const string idFileFolder = @"D:\MassSpecFiles\training\IdResult";
            const string outFileFolder = @"D:\MassSpecFiles\training\FilteredIdResult";

            if (!Directory.Exists(idFileFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
            }


            for (var d = 0; d < TrainSetFileLists.Length; d++)
            {
                var dataset = TrainSetFileLists[d];

                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var idFile = string.Format(@"{0}\{1}_IcTda.tsv", idFileFolder, dataname);
                var outFileName = string.Format(@"{0}\{1}.trainset.tsv", outFileFolder, Path.GetFileNameWithoutExtension(dataset));

                if (File.Exists(outFileName)) continue;

                Console.WriteLine(dataset);

                if (!File.Exists(idFile))
                {
                    idFile = string.Format(@"{0}\{1}_msgfdb_syn.txt", idFileFolder, dataname);

                    if (!File.Exists(idFile))
                    {
                        Console.WriteLine(@"Skipping file since not found: " + idFile);
                        continue;
                    }
                }

                Console.WriteLine(idFile);

                var targetSets = LcMsFeatureTrain.CollectTrainSet(dataset, idFile);
                Console.WriteLine(targetSets.Count);
                

                var writer =
                    new StreamWriter(outFileName);
                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMass\tSequence\tModifications");

                foreach (var prsmSet in targetSets)
                {
                    writer.Write(prsmSet.MinScanNum);
                    writer.Write("\t");
                    writer.Write(prsmSet.MaxScanNum);
                    writer.Write("\t");
                    writer.Write(prsmSet.MinCharge);
                    writer.Write("\t");
                    writer.Write(prsmSet.MaxCharge);
                    writer.Write("\t");
                    writer.Write(prsmSet.Mass);
                    writer.Write("\t");
                    writer.Write(prsmSet[0].Sequence);
                    writer.Write("\t");
                    writer.Write(prsmSet[0].Modifications);
                    //writer.Write("\t")
                    //writer.Write(string.Join(";", prsmSet.Select(prsm => prsm.ScanNum)));
                    writer.Write("\n");
                }

                writer.Close();
            }
        }

        [Test]
        public void ExtractLcMsFeaturesForTrainingSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\FilteredIdResult";
            if (!Directory.Exists(idFileFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
            }
            
            var tolerance = new Tolerance(10);
            var tolerance2 = new Tolerance(20);
            var id = 1;
            
            
            for(var d = 0; d < TrainSetFileLists.Length; d++)
            {
                var dataset = TrainSetFileLists[d];
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var filtedIdResultFile = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var featureResult = string.Format(@"{0}\{1}.ms1ft", idFileFolder, Path.GetFileNameWithoutExtension(dataset));

                if (!File.Exists(dataset))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", dataset);
                    continue;
                }
                if (!File.Exists(filtedIdResultFile))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", filtedIdResultFile);
                    continue;
                }


                var run = PbfLcMsRun.GetLcMsRun(dataset);
                
                
                var targetStatWriter = new StreamWriter(string.Format(@"D:\MassSpecFiles\training\statistics\{0}.tsv", Path.GetFileNameWithoutExtension(dataset)));
                var decoyStatWriter = new StreamWriter(string.Format(@"D:\MassSpecFiles\training\statistics\{0}_decoy.tsv", Path.GetFileNameWithoutExtension(dataset)));
                var writer = new StreamWriter(featureResult);

                writer.Write("Ms2MinScan\tMs2MaxScan\tMs2MinCharge\tMs2MaxCharge\tMs2Mass\t");
                writer.Write("Mass\tMinScan\tMaxScan\tMinCharge\tMaxCharge\tMinTime\tMaxTime\tElution\tGood\n");
                var tsvParser = new TsvFileParser(filtedIdResultFile);

                var featureFinder = new LcMsPeakMatrix(run);
                
                
                for (var i = 0; i < tsvParser.NumData; i++)
                {
                    var minScan = int.Parse(tsvParser.GetData("MinScan")[i]);
                    var maxScan = int.Parse(tsvParser.GetData("MaxScan")[i]);
                    var minCharge = int.Parse(tsvParser.GetData("MinCharge")[i]);
                    var maxCharge = int.Parse(tsvParser.GetData("MaxCharge")[i]);
                    var mass = double.Parse(tsvParser.GetData("Mass")[i]);

                    writer.Write(minScan);
                    writer.Write("\t");
                    writer.Write(maxScan);
                    writer.Write("\t");
                    writer.Write(minCharge);
                    writer.Write("\t");
                    writer.Write(maxCharge);
                    writer.Write("\t");
                    writer.Write(mass);
                    writer.Write("\t");

                    var binNum = featureFinder.Comparer.GetBinNumber(mass);

                    var binMass = featureFinder.Comparer.GetMzAverage(binNum);

                    var binNumList = (mass < binMass) ? new int[] { binNum, binNum - 1, binNum + 1 } : new int[] { binNum, binNum + 1, binNum - 1 };
                    LcMsPeakCluster refinedFeature = null;

                    foreach (var bi in binNumList)
                    {
                        var tempList = new List<LcMsPeakCluster>();
                        var features = featureFinder.FindFeatures(bi);
                        var massTh = (mass < 2000) ? tolerance2.GetToleranceAsTh(mass) : tolerance.GetToleranceAsTh(mass);
                        foreach (var feature in features)
                        {
                            if (Math.Abs(mass - feature.Mass) < massTh) tempList.Add(feature);
                        }

                        //var nHits = 0;
                        var highestAbu = 0d;
                        //var scans = Enumerable.Range(minScan, maxScan - minScan + 1);
                        foreach (var feature in tempList)
                        {
                            //var scans2 = Enumerable.Range(feature.MinScanNum, feature.MaxScanNum - feature.MinScanNum + 1);
                            //var hitScans = scans.Intersect(scans2).Count();
                            if (feature.MinScanNum < 0.5*(minScan + maxScan) &&
                                0.5*(minScan + maxScan) < feature.MaxScanNum)
                            {
                                if (feature.Abundance > highestAbu)
                                {
                                    refinedFeature = feature;
                                    highestAbu = feature.Abundance;
                                }
                            }
                            /*if (hitScans > 0)
                            {
                                refinedFeature = feature;
                                nHits = hitScans;
                            }*/
                        }

                        if (refinedFeature != null) break;
                    }

                    if (refinedFeature != null)
                    {
                        writer.Write(refinedFeature.Mass);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MinScanNum);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MaxScanNum);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MinCharge);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MaxCharge);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MinElutionTime);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MaxElutionTime);
                        writer.Write("\t");
                        writer.Write(refinedFeature.MaxElutionTime - refinedFeature.MinElutionTime);
                        writer.Write("\t");

                        var good = (refinedFeature.MinScanNum <= minScan && refinedFeature.MaxScanNum >= maxScan);
                        writer.Write(good ? 1 : 0);
                        writer.Write("\n");
                        //writer.Write(0); writer.Write("\t");
                        //writer.Write(0); writer.Write("\n");

                        OutputEnvelopPeakStat(id, refinedFeature, targetStatWriter);

                        var chargeRange = featureFinder.GetDetectableMinMaxCharge(refinedFeature.RepresentativeMass, run.MinMs1Mz, run.MaxMs1Mz);
                        refinedFeature.UpdateWithDecoyScore(featureFinder.Ms1Spectra, chargeRange.Item1, chargeRange.Item2);
                        OutputEnvelopPeakStat(id, refinedFeature, decoyStatWriter);
                        id++;
                    }
                    else
                    {
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\t");
                        writer.Write(0);
                        writer.Write("\n");

                    }
                    //var feature = featureFinder.FindLcMsPeakCluster(mass, (int) scan, (int) charge);
                }
                writer.Close();
                targetStatWriter.Close();
                decoyStatWriter.Close();
                Console.WriteLine(dataname);
            }
        }

        private void OutputEnvelopPeakStat(int id, LcMsPeakCluster feature, StreamWriter writer)
        {
            /*
            public double[] EnvelopeDistanceScoreAcrossCharge { get; internal set; }
            public double[] EnvelopeCorrelationScoreAcrossCharge { get; internal set; }
            public double[] EnvelopeIntensityScoreAcrossCharge { get; internal set; }
            public double[] AbundanceDistributionAcrossCharge { get; internal set; }
            public double[] BestCorrelationScoreAcrossCharge { get; private set; }
            public double[] BestDistanceScoreAcrossCharge { get; private set; }
            public double[] BestIntensityScoreAcrossCharge { get; private set; }            
            */

            //for(var charge = feature.MinCharge; charge <= feature.MaxCharge; charge++)
            for (var i = 0; i < 2; i++)
            {
                writer.Write(id);
                writer.Write("\t");

                writer.Write(feature.Mass);
                writer.Write("\t");

                writer.Write(feature.BestCharge[i]);
                writer.Write("\t");

                writer.Write(feature.EnvelopeDistanceScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.EnvelopeCorrelationScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.EnvelopeIntensityScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.BestDistanceScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.BestCorrelationScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.BestIntensityScoreAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.AbundanceDistributionAcrossCharge[i]);
                writer.Write("\t");

                writer.Write(feature.XicCorrelationBetweenBestCharges[0]);
                writer.Write("\t");

                writer.Write(feature.XicCorrelationBetweenBestCharges[1]);
                //writer.Write("\t");

                writer.Write("\n");
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
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, TestRawFile);
            }

            const string TestResultFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01_IcTda.tsv";
            if (!File.Exists(TestResultFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, TestResultFile);
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
                var score = featureFinder.GetMs1EvidenceScore(scan, mass, charge);
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
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, resultFilePath);
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
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const string outFolderPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Output";
            if (!Directory.Exists(outFolderPath))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, outFolderPath);
            }

            //const string specFilePath = @"D:\MassSpecFiles\test\QC_Shew_Intact_4_01Jan15_Bane_C2-14-08-02RZ.raw";
            
            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 5000;
            const int maxThreads = 10;

            var param = new LcMsFeatureFinderInputParameter()
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
            var featureFinder = new LcMsFeatureFinderLauncher(param);
            featureFinder.Run();
        }

        [Test]
        [STAThread]
        public void TestFeatureMapGeneration()
        {
            Console.WriteLine("Testing Working");
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            const string rawFile = @"\\protoapps\UserData\Jungkap\Joshua\testData\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string testFile = @"\\protoapps\UserData\Jungkap\Joshua\FeatureMap\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            const string outputFile = @"D:\MassSpecFiles\training\raw\";

            var map = new LcMsFeatureMap(PbfLcMsRun.GetLcMsRun(rawFile), testFile, 2000, 50000);
            map.SaveImage(outputFile + "test.png");
        }

        [Test]
        public void TestProMexFilter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            
            const string ms1FtPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            if (!File.Exists(ms1FtPath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, ms1FtPath);
            }

            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }



        public static string[] TrainSetFileLists = new string[]
            {
            
                @"D:\MassSpecFiles\bottom_up\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw",  // bottom-up datasets
                @"D:\MassSpecFiles\bottom_up\QC_Shew_13_05_2_17Apr14_Samwise_13-07-17.raw",                

                @"D:\MassSpecFiles\IMER\Dey_IMERblast_01_08May14_Alder_14-01-33.pbf",   
                @"D:\MassSpecFiles\IMER\Dey_IMERblast_02_08May14_Alder_14-01-33.pbf",
                @"D:\MassSpecFiles\IMER\Dey_IMERblast_03_09May14_Alder_14-01-33.pbf",
                
                @"D:\MassSpecFiles\IMER\Diabetes_iPSC_Beta_1_IMER_13May14_Alder_14-01-33.pbf",
                @"D:\MassSpecFiles\IMER\Diabetes_iPSC_Beta_3_IMER_14May14_Alder_14-01-33.pbf",
                
                @"\\protoapps\UserData\Jungkap\peptidome\CPTAC_Peptidome_Test1_P2_13Jan12_Polaroid_11-10-14.pbf", // peptidomics datasets
                @"\\protoapps\UserData\Jungkap\peptidome\CPTAC_Peptidome_Test2_P5-18_13Jan12_Polaroid_11-10-14.pbf",
                @"\\protoapps\UserData\Jungkap\peptidome\CPTAC_Peptidome_Test2_P6-5_13Jan12_Polaroid_11-10-14.pbf",

                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",// top-down datasets
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf",
                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_03_15Jan15_Bane_C2-14-08-02RZ.pbf",

                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_03_15Jan15_Bane_C2-14-08-02RZ.pbf",

                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_3_2Feb15_Bane_C2Column4.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
                @"D:\MassSpecFiles\SBEP\SBEP_STM_001_02222012_Aragon.pbf",
                
                @"D:\MassSpecFiles\test\NewQC_LongSep_29Sep14_141001104925.pbf",
                @"D:\MassSpecFiles\training\raw\YS_Shew_testHCD_CID.pbf",
                //@"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2014_2\yufeng_column_test2.pbf"
            };

    }
    


}
