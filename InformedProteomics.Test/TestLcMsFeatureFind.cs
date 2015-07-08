using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
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
        public void TestIsotopeCount()
        {
            for (var m = 1000; m < 30000; m += 1000)
            {
                var TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(m, 30);
                var n = 2; //Math.Ceiling(TheoreticalEnvelope.Size*0.5);

                var j = 0;
                for (var i = 0; i < TheoreticalEnvelope.IndexOrderByRanking.Length; i++)
                {
                    var k = TheoreticalEnvelope.IndexOrderByRanking[i];
                    if (TheoreticalEnvelope.Isotopes[k].Ratio < 0.3) break;
                    j++;
                }
                Console.WriteLine("{0}\t{1}\t{2}", m, j, TheoreticalEnvelope.Size);

            }
        }
        
        [Test]
        public void TestLcMsFeatureFinder()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            // var outTsvFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, "ms1ft");

            var scoreDataPath = @"D:\MassSpecFiles\training";
            var scorer = new LcMsFeatureLikelihood(scoreDataPath);

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine(@"Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var featureFinder = new LcMsPeakMatrix(run, scorer);
            Console.WriteLine(@"Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

            var container = new LcMsFeatureContainer(featureFinder.Ms1Spectra, scorer);
            var minSearchMassBin = featureFinder.Comparer.GetBinNumber(3000);
            var maxSearchMassBin = featureFinder.Comparer.GetBinNumber(3010);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine(@"Start MS1 feature extraction.");

            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                var clusters = featureFinder.FindFeatures(binNum);
                container.Add(clusters);

                if (binNum > minSearchMassBin && (binNum - minSearchMassBin) % 1000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - minSearchMassBin;
                    var processedPercentage = ((double)processedBins / totalMassBin) * 100;
                    Console.WriteLine(@"Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}", processedPercentage, featureFinder.Comparer.GetMzEnd(binNum), elapsed, container.NumberOfFeatures);
                }
            }

            Console.WriteLine(@"Complete MS1 feature extraction.");
            Console.WriteLine(@" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(@" - Number of extracted features = {0}", container.NumberOfFeatures);

            // write result files
            Console.WriteLine(@"Start selecting mutually independent features from feature network graph");
            var connectedFeatures = container.GetAllConnectedFeatures();
          
            stopwatch.Stop();

            // Start to quantify accurate abundance
            stopwatch.Restart();
            //var quantAnalyzer = new TargetMs1FeatureMatrix(run);
            //var oriResult = new List<Ms1FeatureCluster>();
            //var quantResult = new List<Ms1Feature>();
    

            var featureId = 0;
            var ms1ScanNums = run.GetMs1ScanVector();
            //tsvWriter.WriteLine(GetHeaderString() + "\tQMinScanNum\tQMaxScanNum\tQMinCharge\tQMaxCharge\tQAbundance");
            
            //foreach (var cluster in container.GetFilteredFeatures(connectedFeatures))
            foreach (var feature in container.GetFilteredFeatures(connectedFeatures))
            {
                Console.Write(featureId); Console.Write("\t");
                Console.Write(feature.Mass); Console.Write("\t");
                Console.Write(feature.MinScanNum); Console.Write("\t");
                Console.Write(feature.MaxScanNum); Console.Write("\t");
                Console.Write(feature.MinCharge); Console.Write("\t");
                Console.Write(feature.MaxCharge); Console.Write("\t");
                
                Console.Write(feature.RepresentativeScanNum); Console.Write("\t");
                Console.Write(feature.RepresentativeMz); Console.Write("\t");
                Console.Write(feature.RepresentativeCharge); Console.Write("\t");

                //Console.Write(feature.BestSummedEnvelopeDistance); Console.Write("\t");
                //Console.Write(feature.BestEnvelopeDistance); Console.Write("\t");
                Console.Write(feature.BestDistanceScoreAcrossCharge[0]); Console.Write("\t");
                Console.Write(feature.BestDistanceScoreAcrossCharge[1]); Console.Write("\t");

                Console.Write(feature.BestCorrelationScoreAcrossCharge[0]); Console.Write("\t");
                Console.Write(feature.BestCorrelationScoreAcrossCharge[1]); Console.Write("\t");

                Console.Write(feature.BestIntensityScoreAcrossCharge[0]); Console.Write("\t");
                Console.Write(feature.BestIntensityScoreAcrossCharge[1]); Console.Write("\t");

                Console.Write(feature.AbundanceDistributionAcrossCharge[0]); Console.Write("\t");
                Console.Write(feature.AbundanceDistributionAcrossCharge[1]); Console.Write("\t");

                Console.Write(feature.XicCorrelationBetweenBestCharges[0]); Console.Write("\t");
                Console.Write(feature.XicCorrelationBetweenBestCharges[1]); Console.Write("\t");

                Console.Write(feature.Score); Console.Write("\n");
                featureId++;

            }
        }
        
        
        
        [Test]
        public void CollectTrainingSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\IcTda";
            const string outFileFolder = @"D:\MassSpecFiles\training\refined_set";

            if (!Directory.Exists(idFileFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
                return;
            }

            var rawFileLists = new string[]
            {
                @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf",                
                @"D:\MassSpecFiles\SBEP\SBEP_STM_001_02222012_Aragon.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
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

                if (!File.Exists(idFile))
                {
                    Console.WriteLine(@"Skipping file since not found: " + idFile);
                    continue;
                }

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
            for(var i = 0; i < 2; i++)
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
        public void TestLcMsFeatureLikelihood()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string folder = @"D:\MassSpecFiles\training";
            if (!Directory.Exists(folder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, folder);
                return;
            }

            var likelihood = new LcMsFeatureLikelihood(folder);
            const string idFileFolder = @"D:\MassSpecFiles\training\refined_set";
            if (!Directory.Exists(idFileFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
                return;
            }

            var rawFileLists = new string[]
            {
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
            };
            
            var id = 1;
            foreach (var dataset in rawFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var featureList = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var featureResult = string.Format(@"{0}\{1}.score.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));

                if (!File.Exists(featureList))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", featureList);
                    continue;
                }

                if (!File.Exists(featureResult))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", featureResult);
                    continue;
                }

                var writer = new StreamWriter(featureResult);

                writer.WriteLine("MinScan\tMaxScan\tMinCharge\tMaxCharge\tMinTime\tMaxTime\tElution\tGood");

                var run = PbfLcMsRun.GetLcMsRun(dataset);
                var featureFinder = new LcMsPeakMatrix(run);
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
                    //var scan = 0.5 * (minScan + maxScan);

                    var feature = featureFinder.GetLcMsPeakCluster(mass, (int)charge, minScan, maxScan);
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
                    writer.Write("\t");
                    var score = likelihood.GetScore(feature);
                    writer.Write(score);
                    writer.Write("\n");

                    if (feature.EnvelopeDistanceScoreAcrossCharge == null) continue;

                }
                writer.Close();
                
                Console.WriteLine(dataname);
            }
        }
        /*
        [Test]
        public void TestFeatureElutionProfile()
        {
            const string idFileFolder = @"D:\MassSpecFiles\training\refined_set";
            var rawFileLists = new string[]
            {
                @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf"
            };

            var id = 1;

            foreach (var dataset in rawFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var featureList = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                
                var featureResult = string.Format(@"{0}\{1}.xic.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var writer = new StreamWriter(featureResult);
                var run = PbfLcMsRun.GetLcMsRun(dataset);

                var ms1ScanNums = run.GetMs1ScanVector();
                foreach (var scan in ms1ScanNums)
                {
                    writer.Write(run.GetElutionTime(scan));
                    writer.Write("\t");
                }
                writer.Write("\n");


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

                    var feature = featureFinder.GetLcMsPeakCluster(mass, (int)charge, minScan, maxScan);

                    if (feature == null)
                    {
                        foreach (var x in ms1ScanNums)
                        {
                            writer.Write(0);
                            writer.Write("\t");
                        }
                        writer.Write("\n");
                        continue;
                    }

                    var xic = featureFinder.GetEntireXic(feature);
                    foreach (var x in xic)
                    {
                        writer.Write(x);
                        writer.Write("\t");                        
                    }
                    writer.Write("\n");

                }
                writer.Close();
                Console.WriteLine(dataname);
            }
        }
        */

        [Test]
        public void TestFeatureFindingWithMassBinning()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\refined_set";
            if (!Directory.Exists(idFileFolder))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
                return;
            }

            var rawFileLists = new string[]
            {
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",                                
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_3_2Feb15_Bane_C2Column4.pbf",                
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",                
                @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf",                
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"D:\MassSpecFiles\SBEP\SBEP_STM_001_02222012_Aragon.pbf",
                @"D:\MassSpecFiles\training\raw\YS_Shew_testHCD_CID.pbf",
                @"D:\MassSpecFiles\test\NewQC_LongSep_29Sep14_141001104925.pbf",
            };

            var tolerance = new Tolerance(10);
            var id = 1;
            foreach (var dataset in rawFileLists)
            {
                if (!File.Exists(dataset))
                {
                    Console.WriteLine(@"Warning: Skipping since file not found: {0}", dataset);
                    continue;
                }

                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var featureList = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder,
                    Path.GetFileNameWithoutExtension(dataset));
                
                var run = PbfLcMsRun.GetLcMsRun(dataset);
                var featureFinder = new LcMsPeakMatrix(run);
                var tsvParser = new TsvFileParser(featureList);
                var featureResult = string.Format(@"{0}\{1}.ms1ft", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var targetStatWriter = new StreamWriter(string.Format(@"D:\MassSpecFiles\training\stats\{0}.tsv", Path.GetFileNameWithoutExtension(dataset)));
                var decoyStatWriter = new StreamWriter(string.Format(@"D:\MassSpecFiles\training\stats\{0}_decoy.tsv", Path.GetFileNameWithoutExtension(dataset)));
                var writer = new StreamWriter(featureResult);

                writer.Write("Ms2MinScan\tMs2MaxScan\tMs2MinCharge\tMs2MaxCharge\tMs2Mass\t");
                writer.Write("Mass\tMinScan\tMaxScan\tMinCharge\tMaxCharge\tMinTime\tMaxTime\tElution\tGood\tCorr\tDist");
                writer.WriteLine("\tseedDist\tseedCorr\tseedPoisson\tseedRanksum");

                for (var i = 0; i < tsvParser.NumData; i++)
                {
                    var minScan = int.Parse(tsvParser.GetData("MinScan")[i]);
                    var maxScan = int.Parse(tsvParser.GetData("MaxScan")[i]);
                    var minCharge = int.Parse(tsvParser.GetData("MinCharge")[i]);
                    var maxCharge = int.Parse(tsvParser.GetData("MaxCharge")[i]);
                    var mass = double.Parse(tsvParser.GetData("Mass")[i]);
                    var sequence = tsvParser.GetData("Sequence")[i];
                    var composition = tsvParser.GetData("Composition")[i];
                    var charge = 0.5*(minCharge + maxCharge);
                    var scan = 0.5*(minScan + maxScan);

                    writer.Write(minScan); writer.Write("\t");
                    writer.Write(maxScan); writer.Write("\t");
                    writer.Write(minCharge); writer.Write("\t");
                    writer.Write(maxCharge); writer.Write("\t");
                    writer.Write(mass); writer.Write("\t");

                    var binNum = featureFinder.Comparer.GetBinNumber(mass);
                    
                    var binNumList = new int[] {binNum, binNum - 1, binNum + 1};
                    LcMsPeakCluster refinedFeature = null;

                    foreach(var bi in binNumList)
                    {
                        var tempList = new List<LcMsPeakCluster>();
                        var features = featureFinder.FindFeatures(bi);
                        var massTh = tolerance.GetToleranceAsTh(mass);
                        foreach (var feature in features)
                        {
                            if (Math.Abs(mass - feature.Mass) < massTh) tempList.Add(feature);
                        }
                        
                        var nHits = 0;
                        var scans = Enumerable.Range(minScan, maxScan - minScan + 1);
                        foreach (var feature in tempList)
                        {
                            var scans2 = Enumerable.Range(feature.MinScanNum, feature.MaxScanNum - feature.MinScanNum + 1);
                            var hitScans = scans.Intersect(scans2).Count();
                            if (hitScans > nHits)
                            {
                                refinedFeature = feature;
                                nHits = hitScans;
                            }
                        }

                        if (nHits > 0) break;
                    }
                    
                    if (refinedFeature != null)
                    {
                        writer.Write(refinedFeature.Mass); writer.Write("\t");
                        writer.Write(refinedFeature.MinScanNum); writer.Write("\t");
                        writer.Write(refinedFeature.MaxScanNum); writer.Write("\t");
                        writer.Write(refinedFeature.MinCharge); writer.Write("\t");
                        writer.Write(refinedFeature.MaxCharge); writer.Write("\t");
                        writer.Write(refinedFeature.MinElutionTime); writer.Write("\t");
                        writer.Write(refinedFeature.MaxElutionTime); writer.Write("\t");
                        writer.Write(refinedFeature.MaxElutionTime - refinedFeature.MinElutionTime); writer.Write("\t");

                        var good = (refinedFeature.MinScanNum <= minScan && refinedFeature.MaxScanNum >= maxScan);

                        writer.Write(good ? 1 : 0); writer.Write("\n");
                        //writer.Write(0); writer.Write("\t");
                        //writer.Write(0); writer.Write("\n");
                        
                        OutputEnvelopPeakStat(id, refinedFeature, targetStatWriter);
                        
                        refinedFeature.UpdateWithDecoyScore(featureFinder.Ms1Spectra);
                        OutputEnvelopPeakStat(id, refinedFeature, decoyStatWriter);
                        id++;

                    }
                    else
                    {
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\t");
                        writer.Write(0); writer.Write("\n");

                    }
                    //var feature = featureFinder.FindLcMsPeakCluster(mass, (int) scan, (int) charge);
                }
                writer.Close();
                targetStatWriter.Close();
                decoyStatWriter.Close();
                Console.WriteLine(dataname);
            }
            
        }
    }
}
