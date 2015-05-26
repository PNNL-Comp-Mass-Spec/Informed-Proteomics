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
    public class TestProMex
    {

        [Test]
        public void CollectTrainingSet()
        {
            const string idFileFolder = @"\\protoapps\UserData\Jungkap\TrainingSet";
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
            const string rawFile = @"D:\Test\Quant\spike_in\raw\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ.pbf";
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
            for (var m = 500; m < 5000; m += 100)
            {
                var theoreticalEnvelope = new IsotopeList(m, 30, 0.1);    

                Console.WriteLine("{0}\t{1}", m, theoreticalEnvelope.Count);
            }
            
        }

        [Test]
        public void TestMs1EvidenceScore()
        {
            const string TestRawFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01.pbf";
            const string TestResultFile = @"\\protoapps\UserData\Jungkap\Lewy\Lewy_intact_01_IcTda.tsv";

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
            var featureDir = @"D:\MassSpecFiles\Plasma";
            var rawDir = @"\\protoapps\UserData\Jungkap\plasma\raw";
            var outFile = @"D:\MassSpecFiles\Plasma\aligned_features.tsv";
            var dataset = new string[4]
            {
                "UC4_Intact_plasmaTest_41_6May15_Bane_14-09-01RZ",
                "UC4_Intact_plasmaTest_144_6May15_Bane_14-09-01RZ",
                "UC4_Intact_plasmaTest_90_6May15_Bane_14-09-01RZ",
                "UC4_Intact_plasmaTest_92_6May15_Bane_14-09-01RZ"
            };

            /*
            var featureDir = @"D:\Test\Quant\spike_in\ms1ft";
            var rawDir = @"\\protoapps\UserData\Jungkap\Quant";
            var outFile = @"D:\Test\Quant\spike_in\aligned_features.tsv";

            var dataset = new string[15]
            {
                "CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ",
                
                "CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ",
                "CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ",
                
                "CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ",
                "CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ",
            };
            */

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
            /*
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
            */
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
            //align.AlignNetAll(maxDatasetIndex);

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

                var features = new Ms1Feature[dataset.Length];
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
            /*
            for (var i = 0; i < nDataset; i++)
            {
                var wt = new StreamWriter(string.Format(@"{0}\{1}.fid2gid", featureDir, dataset[i]));
                for(var j = 0; j < fid2gid[i].Length; j++) wt.WriteLine(fid2gid[i][j]);
                wt.Close();
            }
            */
            
        }

        
        private void OutputEnvelopPeakStat(Ms1Feature feature, int baseCharge, LcMsRun run, IsotopeList isotopes, StreamWriter writer)
        {
            var envelopes = feature.EnvelopeList;
            var ms1ScanNums = run.GetMs1ScanVector();
            var minEt = feature.MinElutionTime;
            var maxEt = feature.MaxElutionTime;
            var etLen = maxEt - minEt;
            
            foreach (var env in envelopes)
            {
                var charge = env.Row + baseCharge;
                var et = run.GetElutionTime(ms1ScanNums[env.Col]);
                var net = (et - minEt)/etLen;

                for(var k = 0; k < isotopes.Count; k++)
                {
                    var peak = env.Peaks[k];
                    if (peak == null) continue;

                    var isoRanking = isotopes.IsotopeRanking[k];
                    var isoIdx = isotopes.SortedIndexByIntensity[k];
                    writer.Write("{0:0.0}\t{1}\t{2}\t{3}\t", feature.Mass, charge, net, isoRanking);

                    //if (feature.Mass < 15000) writer.Write(peak.LocalRankingForLowMass);
                    //else writer.Write(peak.LocalRankingForHighMass);

                    writer.Write("\n");
                }
            }

        }

        [Test]
        public void FeatureAnalysisInTrainingSet()
        {
            const string idFileFolder = @"D:\MassSpecFiles\training\IcTda";
            var trainingSetList = File.ReadAllLines(@"D:\MassSpecFiles\training\training_datasets.txt");
            var trainSet = trainingSetList.Where(set => !set.StartsWith("#")).ToList();

            foreach (var dataset in trainSet)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var featureList = string.Format(@"D:\MassSpecFiles\training\target\{0}.target.tsv",
                    Path.GetFileNameWithoutExtension(dataset));

                var featureResult = string.Format(@"D:\MassSpecFiles\training\target\{0}.result.tsv",
                    Path.GetFileNameWithoutExtension(dataset));

                var featureStat = string.Format(@"D:\MassSpecFiles\training\stats\{0}.tsv",
                                    Path.GetFileNameWithoutExtension(dataset));

                var statWriter = new StreamWriter(featureStat);
                var writer = new StreamWriter(featureResult);

                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMinTime\tMaxTime\tElution\tGood");

                var run = PbfLcMsRun.GetLcMsRun(dataset);
                var featureFinder = new TargetMs1FeatureMatrix(run);
                var tsvParser = new TsvFileParser(featureList);

                for (var i = 0; i < tsvParser.NumData; i++)
                {
                    //var featureDist = string.Format(@"D:\MassSpecFiles\training\{0}\{1}.txt", Path.GetFileNameWithoutExtension(dataset), i);
                    //var distWriter = new StreamWriter(featureDist);

                    var minScan = int.Parse(tsvParser.GetData("MinScan")[i]);
                    var maxScan = int.Parse(tsvParser.GetData("MaxScan")[i]);
                    var minCharge = int.Parse(tsvParser.GetData("MinCharge")[i]);
                    var maxCharge = int.Parse(tsvParser.GetData("MaxCharge")[i]);
                    var mass = double.Parse(tsvParser.GetData("Mass")[i]);
                    var sequence = tsvParser.GetData("Sequence")[i];
                    var composition = tsvParser.GetData("Composition")[i];

                    var charge = 0.5*(minCharge + maxCharge);
                    var scan = 0.5*(minScan + maxScan);

                    var target = new TargetFeature(mass, (int) charge, (int) scan)
                    {
                        MinCharge = minCharge,
                        MaxCharge = maxCharge,
                        MinScanNum = minScan,
                        MaxScanNum = maxScan,
                        Sequence = sequence,
                        Composition = composition
                    };
                    var ms1Feature = featureFinder.FindMs1Feature(target);

                    writer.Write(ms1Feature.MinScanNum);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MaxScanNum);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MinCharge);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MaxCharge);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MinElutionTime);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MaxElutionTime);
                    writer.Write("\t");
                    writer.Write(ms1Feature.MaxElutionTime - ms1Feature.MinElutionTime);
                    writer.Write("\t");

                    var good = (ms1Feature.MinScanNum <= target.MinScanNum && ms1Feature.MaxScanNum >= target.MaxScanNum);
                    //&& ms1Feature.MinCharge <= target.MinCharge && ms1Feature.MaxCharge >= target.MaxCharge);

                    writer.Write(good ? 1 : 0);
                    writer.Write("\n");

                    OutputEnvelopPeakStat(ms1Feature, featureFinder.MinScanCharge, run, featureFinder.TheoreticalEnvelope, statWriter);
                    //Console.WriteLine("{0}\t{1}", i, composition);
                    //distWriter.Write(ms1Feature.Desc);
                    //distWriter.Close();
                }
                writer.Close();
                statWriter.Close();
                Console.WriteLine(dataname);
            }
        }

        [Test]
        public void ProMexDebug()
        {
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
