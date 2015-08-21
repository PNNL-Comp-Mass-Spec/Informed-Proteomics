﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
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


            var threads = 10;

            var ms2Scans = Enumerable.Range(1, 7892).ToList();
            var nScansPerThread = (int)ms2Scans.Count/threads;



            Math.Ceiling((double) ms2Scans.Count/(double) nScansPerThread);

            for (var m = 500; m < 3000; m += 100)
            {
                var TheoreticalEnvelope = new TheoreticalIsotopeEnvelope(m, 30);
                //var n = 2; //Math.Ceiling(TheoreticalEnvelope.Size*0.5);

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
        public void TestLcMsFeatureXic()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFile = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf";
            //const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var scorer = new LcMsFeatureLikelihood();
            var featureFinder = new LcMsPeakMatrix(run, scorer);
            var feature = featureFinder.GetLcMsPeakCluster(2388.278, 4, 3774, 3907);

            //feature = featureFinder.GetLcMsPeakCluster(8151.3706, 7, 13, 4201, 4266);

            //feature = featureFinder.GetLcMsPeakCluster(8151.41789, 7, 13, 2861, 2941);
            
            var ms1ScanToIndex = run.GetMs1ScanNumToIndex();
            var minCol = ms1ScanToIndex[feature.MinScanNum];
            var maxCol = ms1ScanToIndex[feature.MaxScanNum];
            //var minRow = feature.MinCharge - LcMsPeakMatrix.MinScanCharge;
            //var maxRow = feature.MaxCharge - LcMsPeakMatrix.MinScanCharge;
            
            Console.WriteLine("---------------------------------------------------------------");
            for (var i = 0; i < feature.Envelopes.Length; i++)
            {
                for (var j = 0; j < feature.Envelopes[i].Length; j++)
                {
                    Console.Write(feature.Envelopes[i][j] != null ? feature.Envelopes[i][j].PearsonCorrelation : 0);
                    Console.Write("\t");
                }
                Console.Write("\n");
            }
            Console.WriteLine("---------------------------------------------------------------");
            for (var i = 0; i < feature.Envelopes.Length; i++)
            {
                for (var j = 0; j < feature.Envelopes[i].Length; j++)
                {
                    Console.Write(feature.Envelopes[i][j] != null ? feature.Envelopes[i][j].BhattacharyyaDistance : 0);
                    Console.Write("\t");
                }
                Console.Write("\n");
            }

            Console.WriteLine("---------------------------------------------------------------");
            for (var i = 0; i < feature.Envelopes.Length; i++)
            {
                for (var j = 0; j < feature.Envelopes[i].Length; j++)
                {
                    Console.Write(feature.Envelopes[i][j] != null ? feature.Envelopes[i][j].Abundance : 0);
                    Console.Write("\t");
                }
                Console.Write("\n");
            }
            
            



        }

        [Test]
        public void TestLcMsFeatureFinder()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            //const string rawFile = @"D:\MassSpecFiles\CompRef\CPTAC_Intact_CR_Pool_2_25Jun15_Bane_15-02-02RZ.pbf";
            //const string rawFile = @"D:\MassSpecFiles\IMER\Dey_IMERblast_01_08May14_Alder_14-01-33.pbf";
            const string rawFile = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_3\MZ20150729FG_WT1.pbf";

            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            // var outTsvFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, "ms1ft");
            //var scoreDataPath = @"D:\MassSpecFiles\training";
            var scorer = new LcMsFeatureLikelihood();
            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine(@"Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var featureFinder = new LcMsPeakMatrix(run, scorer);
            Console.WriteLine(@"Complete loading MS1 data. Elapsed Time = {0:0.000} sec",
                (stopwatch.ElapsedMilliseconds)/1000.0d);

            var container = new LcMsFeatureContainer(featureFinder.Ms1Spectra, scorer, new LcMsFeatureMergeComparer(new Tolerance(10)));
            var minSearchMassBin = featureFinder.Comparer.GetBinNumber(14138.7);
            var maxSearchMassBin = featureFinder.Comparer.GetBinNumber(14138.7);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine(@"Start MS1 feature extraction.");

            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                var clusters = featureFinder.FindFeatures(binNum);
                container.Add(clusters);

                if (binNum > minSearchMassBin && (binNum - minSearchMassBin)%1000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds)/1000.0d;
                    var processedBins = binNum - minSearchMassBin;
                    var processedPercentage = ((double) processedBins/totalMassBin)*100;
                    Console.WriteLine(
                        @"Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}",
                        processedPercentage, featureFinder.Comparer.GetMzEnd(binNum), elapsed,
                        container.NumberOfFeatures);
                }
            }

            Console.WriteLine(@"Complete MS1 feature extraction.");
            Console.WriteLine(@" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds)/1000.0d);
            Console.WriteLine(@" - Number of extracted features = {0}", container.NumberOfFeatures);

            // write result files
            Console.WriteLine(@"Start selecting mutually independent features from feature network graph");
            

            stopwatch.Stop();

            // Start to quantify accurate abundance
            stopwatch.Restart();
            //var quantAnalyzer = new TargetMs1FeatureMatrix(run);
            //var oriResult = new List<Ms1FeatureCluster>();
            //var quantResult = new List<Ms1Feature>();

            var featureId = 0;
            var ms1ScanNums = run.GetMs1ScanVector();
            //tsvWriter.WriteLine(GetHeaderString() + "\tQMinScanNum\tQMaxScanNum\tQMinCharge\tQMaxCharge\tQAbundance");
            
            var filteredFeatures = container.GetFilteredFeatures(featureFinder);
            foreach (var feature in filteredFeatures)
            {
                Console.Write(featureId);
                Console.Write("\t");
                Console.Write(feature.Mass);
                Console.Write("\t");
                Console.Write(feature.MinScanNum);
                Console.Write("\t");
                Console.Write(feature.MaxScanNum);
                Console.Write("\t");
                Console.Write(feature.MinCharge);
                Console.Write("\t");
                Console.Write(feature.MaxCharge);
                Console.Write("\t");

                Console.Write(feature.RepresentativeScanNum);
                Console.Write("\t");
                Console.Write(feature.RepresentativeMz);
                Console.Write("\t");
                Console.Write(feature.RepresentativeCharge);
                Console.Write("\t");

                //Console.Write(feature.BestSummedEnvelopeDistance); Console.Write("\t");
                //Console.Write(feature.BestEnvelopeDistance); Console.Write("\t");
                Console.Write(feature.BestDistanceScoreAcrossCharge[0]);
                Console.Write("\t");
                Console.Write(feature.BestDistanceScoreAcrossCharge[1]);
                Console.Write("\t");

                Console.Write(feature.BestCorrelationScoreAcrossCharge[0]);
                Console.Write("\t");
                Console.Write(feature.BestCorrelationScoreAcrossCharge[1]);
                Console.Write("\t");

                Console.Write(feature.BestIntensityScoreAcrossCharge[0]);
                Console.Write("\t");
                Console.Write(feature.BestIntensityScoreAcrossCharge[1]);
                Console.Write("\t");

                Console.Write(feature.AbundanceDistributionAcrossCharge[0]);
                Console.Write("\t");
                Console.Write(feature.AbundanceDistributionAcrossCharge[1]);
                Console.Write("\t");

                Console.Write(feature.XicCorrelationBetweenBestCharges[0]);
                Console.Write("\t");
                Console.Write(feature.XicCorrelationBetweenBestCharges[1]);
                Console.Write("\t");

                Console.Write(feature.Score);
                Console.Write("\n");
                featureId++;

            }
        }



    }
}
