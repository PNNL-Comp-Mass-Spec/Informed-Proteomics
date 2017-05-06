using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Clustering;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.FeatureFinding.Scoring;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests
{
    [TestFixture]
    public class TestLcMsFeatureFind
    {
        [Test]
        [Category("PNL_Domain")]
        public void TestMaxEntDeconvoluter()
        {
            const string rawFileFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_214\2015_4";
            const string fname = "WHIM2_LoHi_T2DD_HCD_GF07_02";
            string rawFile = string.Format(@"{0}\{1}.pbf", rawFileFolder, fname);
            string ms1ft = string.Format(@"\\protoapps\UserData\Jungkap\CompRef\lowRes\{0}.ms1ft", fname);

            var run = PbfLcMsRun.GetLcMsRun(rawFile, 1.4826, 0);
            var ms1ScanNums = run.GetMs1ScanVector();

            var featureFinder = new LcMsPeakMatrixLowResolution(run);

            foreach (var scan in ms1ScanNums)
            {
                var fts = featureFinder.DetectMs1Features(scan);
                //Console.WriteLine("{0}\t{1}",scan, fts.Count);
            }

            var features = featureFinder.GetLcMsFeatures();

            var writer = new StreamWriter(ms1ft);
            var id = 1;
            writer.WriteLine("FeatureID\tMinScan\tMaxScan\tMinCharge\tMaxCharge\tMonoMass\tAbundance\tRepScan\tMaxElutionTime\tElutionLength\tLikelihoodRatio");
            foreach (var feature in features.OrderBy(f => f.Mass))
            {
                writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}", id, feature.MinScanNum, feature.MaxScanNum, feature.MinCharge,
                    feature.MaxCharge, feature.Mass,  feature.Abundance, feature.RepresentativeScanNum, feature.MinElutionTime, feature.MaxElutionTime, 0);
                id++;
            }
            writer.Close();
        }

        [Test]
        public void TestLcMsFeatureXic()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFile = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf";
            //const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFile);
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
            Utils.ShowStarting(methodName);

            var rawFile = Utils.GetTestFile(methodName, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.pbf");

            // var outTsvFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, "ms1ft");
            //var scoreDataPath = @"D:\MassSpecFiles\training";
            var scorer = new LcMsFeatureLikelihood();
            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine(@"Start loading MS1 data from {0}", rawFile.FullName);

            var run = PbfLcMsRun.GetLcMsRun(rawFile.FullName);
            var featureFinder = new LcMsPeakMatrix(run, scorer);
            Console.WriteLine(@"Complete loading MS1 data. Elapsed Time = {0:0.000} sec",
                (stopwatch.ElapsedMilliseconds)/1000.0d);

            var container = new LcMsFeatureContainer(featureFinder.Ms1Spectra, scorer, new LcMsFeatureMergeComparer(new Tolerance(10)));
            var minSearchMassBin = featureFinder.Comparer.GetBinNumber(11180.33677);
            var maxSearchMassBin = featureFinder.Comparer.GetBinNumber(11180.33677);
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
