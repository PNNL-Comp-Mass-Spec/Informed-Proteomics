using System;
using System.IO;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestLcMsFeatureFind
    {

        [Test]
        public void TestCollectTrainingSet()
        {
            const string idFileFolder = @"D:\MassSpecFiles\training\IcTda";
            const string outFileFolder = @"D:\MassSpecFiles\training\refined_set";
            var rawFileLists = new string[]
            {
                @"D:\MassSpecFiles\SBEP\SBEP_STM_001_02222012_Aragon.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
                @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf",
                @"D:\MassSpecFiles\training\raw\YS_Shew_testHCD_CID.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_3_2Feb15_Bane_C2Column4.pbf",
                @"D:\MassSpecFiles\test\NewQC_LongSep_29Sep14_141001104925.pbf",                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf"
            };
            
            foreach (var dataset in rawFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var idFile = string.Format(@"{0}\{1}_IcTda.tsv", idFileFolder, dataname);

                Console.WriteLine(dataset);
                Console.WriteLine(idFile);
                var targetSets = LcMsFeatureTrain.CollectTrainSet(dataset, idFile, 0.005);
                Console.WriteLine(targetSets.Count);
                var outFileName = string.Format(@"{0}\{1}.trainset.tsv", outFileFolder, Path.GetFileNameWithoutExtension(dataset));
                if (File.Exists(outFileName)) continue;

                var writer =
                    new StreamWriter(outFileName);
                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMass\tSequence\tModifications\tComposition");

                foreach (var feature in targetSets)
                {
                    if (feature.Rows.Count < 2) continue;

                    writer.Write(feature.MinScanNum);
                    writer.Write("\t");
                    writer.Write(feature.MaxScanNum);
                    writer.Write("\t");
                    writer.Write(feature.MinCharge);
                    writer.Write("\t");
                    writer.Write(feature.MaxCharge);
                    writer.Write("\t");
                    writer.Write(feature.Mass);
                    writer.Write("\t");
                    writer.Write(feature.Sequence);
                    writer.Write("\t");
                    writer.Write(feature.Modifications);
                    writer.Write("\t");
                    writer.Write(feature.Composition);
                    writer.Write("\n");
                }

                writer.Close();
            }
        }

        /*
        private void OutputEnvelopPeakStat(LcMsFeature feature, int baseCharge, LcMsRun run, IsotopeList isotopes, StreamWriter writer)
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
                var net = (et - minEt) / etLen;

                for (var k = 0; k < isotopes.Count; k++)
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

        }*/

        [Test]
        public void TestFeatureExtraction()
        {
            var dataset = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\NCR_Int_C1_plus2_30k_21May15_Bane_14-09-01RZ.pbf";
            var run = PbfLcMsRun.GetLcMsRun(dataset);
            var featureFinder = new LcMsPeakMatrix(run, 2, 60, 0);

            var mass = 62947.4149;
            var scan = 4211;
            var charge = 37;

            var feature = featureFinder.FindLcMsPeakCluster(mass, scan, charge);

            

        }

        [Test]
        public void FeatureAnalysisInTrainingSet()
        {
            const string idFileFolder = @"D:\MassSpecFiles\training\refined_set";
            var rawFileLists = new string[]
            {
                @"D:\MassSpecFiles\SBEP\SBEP_STM_001_02222012_Aragon.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
                @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf",
                /*@"D:\MassSpecFiles\training\raw\YS_Shew_testHCD_CID.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_3_2Feb15_Bane_C2Column4.pbf",
                @"D:\MassSpecFiles\test\NewQC_LongSep_29Sep14_141001104925.pbf",                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf"*/
            };


            foreach (var dataset in rawFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var featureList = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var featureResult = string.Format(@"{0}\{1}.feature.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var featureStat = string.Format(@"D:\MassSpecFiles\training\stats\{0}.tsv", Path.GetFileNameWithoutExtension(dataset));

                //if (File.Exists(featureResult)) continue;

                var statWriter = new StreamWriter(featureStat);
                var writer = new StreamWriter(featureResult);

                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMinTime\tMaxTime\tElution\tGood");

                var run = PbfLcMsRun.GetLcMsRun(dataset);
                var featureFinder = new LcMsPeakMatrix(run, 2, 60, 0);
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

                    var charge = 0.5 * (minCharge + maxCharge);
                    var scan = 0.5 * (minScan + maxScan);

                    var feature = featureFinder.FindLcMsPeakCluster(mass, (int) scan, (int) charge);

                    writer.Write(feature.MinScanNum);
                    writer.Write("\t");
                    writer.Write(feature.MaxScanNum);
                    writer.Write("\t");
                    writer.Write(feature.MinCharge);
                    writer.Write("\t");
                    writer.Write(feature.MaxCharge);
                    writer.Write("\t");
                    writer.Write(feature.MinElutionTime);
                    writer.Write("\t");
                    writer.Write(feature.MaxElutionTime);
                    writer.Write("\t");
                    writer.Write(feature.MaxElutionTime - feature.MinElutionTime);
                    writer.Write("\t");

                    var good = (feature.MinScanNum <= minScan && feature.MaxScanNum >= maxScan);

                    writer.Write(good ? 1 : 0);
                    writer.Write("\n");

                    //OutputEnvelopPeakStat(ms1Feature, featureFinder.MinScanCharge, run, feature.Envelopes[0].TheoreticalEnvelope, statWriter);
                    //Console.WriteLine("{0}\t{1}", i, composition);
                    //distWriter.Write(ms1Feature.Desc);
                    //distWriter.Close();
                }
                writer.Close();
                statWriter.Close();
                Console.WriteLine(dataname);
            }
        }




    }
}
