using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.FeatureFinding;
using InformedProteomics.FeatureFinding.Clustering;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.FeatureFinding.Graphics;
using InformedProteomics.FeatureFinding.Scoring;
using InformedProteomics.FeatureFinding.Training;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestProMex
    {
        // Ignore Spelling: msgfdb, prsm, Lewy, Jungkap, Samwise, Dey, peptidomics, IMERblast, peptidome, LongSep, yufeng

        private string mFeatureMapPbfFile;
        private string mFeatureMapResultsFile;

        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        [Category("Local_Testing")]
        public void CollectTrainingSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\IdResult";
            const string outFileFolder = @"D:\MassSpecFiles\training\FilteredIdResult";

            if (!Directory.Exists(idFileFolder))
            {
                Assert.Ignore("Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
            }

            foreach (var dataset in TrainSetFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var idFile = string.Format(@"{0}\{1}_IcTda.tsv", idFileFolder, dataname);
                var outFileName = string.Format(@"{0}\{1}.trainset.tsv", outFileFolder, Path.GetFileNameWithoutExtension(dataset));

                if (File.Exists(outFileName))
                {
                    continue;
                }

                Console.WriteLine(dataset);

                if (!File.Exists(idFile))
                {
                    idFile = string.Format(@"{0}\{1}_msgfdb_syn.txt", idFileFolder, dataname);

                    if (!File.Exists(idFile))
                    {
                        Console.WriteLine("Skipping file since not found: " + idFile);
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
        [Category("Local_Testing")]
        public void ExtractLcMsFeaturesForTrainingSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\FilteredIdResult";
            if (!Directory.Exists(idFileFolder))
            {
                Assert.Ignore("Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
            }

            var tolerance = new Tolerance(10);
            var tolerance2 = new Tolerance(20);
            var id = 1;

            foreach (var dataset in TrainSetFileLists)
            {
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var filtedIdResultFile = string.Format(@"{0}\{1}.trainset.tsv", idFileFolder, Path.GetFileNameWithoutExtension(dataset));
                var featureResult = string.Format(@"{0}\{1}.ms1ft", idFileFolder, Path.GetFileNameWithoutExtension(dataset));

                if (!File.Exists(dataset))
                {
                    Console.WriteLine("Warning: Skipping since file not found: {0}", dataset);
                    continue;
                }
                if (!File.Exists(filtedIdResultFile))
                {
                    Console.WriteLine("Warning: Skipping since file not found: {0}", filtedIdResultFile);
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

                    var binNumList = (mass < binMass) ? new[] { binNum, binNum - 1, binNum + 1 } : new[] { binNum, binNum + 1, binNum - 1 };
                    LcMsPeakCluster refinedFeature = null;

                    foreach (var bi in binNumList)
                    {
                        var tempList = new List<LcMsPeakCluster>();
                        var features = featureFinder.FindFeatures(bi);
                        var massTh = (mass < 2000) ? tolerance2.GetToleranceAsMz(mass) : tolerance.GetToleranceAsMz(mass);
                        foreach (var feature in features)
                        {
                            if (Math.Abs(mass - feature.Mass) < massTh)
                            {
                                tempList.Add(feature);
                            }
                        }

                        //var nHits = 0;
                        var highestAbu = 0d;
                        //var scans = Enumerable.Range(minScan, maxScan - minScan + 1);
                        foreach (var feature in tempList)
                        {
                            //var scans2 = Enumerable.Range(feature.MinScanNum, feature.MaxScanNum - feature.MinScanNum + 1);
                            //var hitScans = scans.Intersect(scans2).Count();
                            if (feature.MinScanNum < 0.5 * (minScan + maxScan) &&
                                0.5 * (minScan + maxScan) < feature.MaxScanNum)
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

                        if (refinedFeature != null)
                        {
                            break;
                        }
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

        private void OutputEnvelopPeakStat(int id, LcMsPeakCluster feature, TextWriter writer)
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
        [Category("PNL_Domain")]
        public void TestMs1EvidenceScore()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var testRawFile = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\Lewy_intact_01.pbf");
            if (!File.Exists(testRawFile))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, testRawFile);
            }

            var testResultFile = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\TestOutput\Lewy_intact_01_IcTda.tsv");
            if (!File.Exists(testResultFile))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, testResultFile);
            }

            var run = PbfLcMsRun.GetLcMsRun(testRawFile);
            var tsvParser = new TsvFileParser(testResultFile);
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
        [TestCase(0, 27)]
        [TestCase(1, 31)]
        [TestCase(2, 30)]
        [TestCase(3, 30)]
        [TestCase(4, 29)]
        [TestCase(5, 29)]
        [TestCase(6, 29)]
        [TestCase(7, 28)]
        [TestCase(8, 28)]
        [TestCase(9, 28)]
        [TestCase(10, 28)]
        [TestCase(11, 28)]
        [TestCase(12, 28)]
        [TestCase(13, 27)]
        [TestCase(14, 27)]
        [TestCase(15, 27)]
        [TestCase(16, 27)]
        [TestCase(17, 27)]
        [TestCase(18, 27)]
        [TestCase(19, 27)]
        [TestCase(20, 27)]
        [TestCase(22, 27)]
        [TestCase(24, 27)]
        [TestCase(25, 26)]
        [TestCase(26, 26)]
        [TestCase(28, 26)]
        [TestCase(30, 26)]
        [TestCase(32, 26)]
        [TestCase(50, 25)]
        [TestCase(75, 25)]
        [TestCase(100, 24)]
        [TestCase(200, 24)]
        public void TestConvertPPMResolutionToBinCount(int ppmResolution, int expectedBinCount)
        {
            var binCount = LcMsFeatureFinderInputParameters.GetBitCountForPPMResolution(ppmResolution);

            Console.WriteLine("Resolution of {0} ppm converts to {1} bins", ppmResolution, expectedBinCount);
            Assert.AreEqual(expectedBinCount, binCount);
        }

        [Test]
        [TestCase(0.001, 70)]
        [TestCase(0.005, 74)]
        [TestCase(0.01, 80)]
        [TestCase(0.02, 83)]
        [TestCase(0.05, 98)]
        public void TestParsePtms(double qValueThreshold, int expectedPTMListCount)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var resultFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt_IcTda.tsv");
            var resultFile = Utils.GetTestFile(methodName, resultFilePath);

            var parser = new TsvFileParser(resultFile.FullName);
            var sequences = parser.GetData("Sequence");
            var modifications = parser.GetData("Modifications");
            var compositions = parser.GetData("Composition");
            var scanNums = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var qValues = parser.GetData("QValue").Select(s => Convert.ToDouble(s)).ToArray();
            var nMacthed = parser.GetData("#MatchedFragments");
            var aaSet = new AminoAcidSet();
            var ptmList = new List<Tuple<int, double, double>>();

            var trimChars = new[] { '"', ' ' };
            var filterPassingResults = 0;

            for (var i = 0; i < parser.NumData; i++)
            {
                if (qValues[i] > qValueThreshold)
                {
                    continue;
                }

                filterPassingResults++;

                var seq = new Sequence(sequences[i], aaSet);
                var sequenceComp = seq.Composition + Composition.H2O;

                var modComposition = Composition.Zero;
                var modsStr = modifications[i].Trim(trimChars);
                if (modsStr.Length == 0)
                {
                    continue;
                }

                var mods = modsStr.Split(',');
                foreach (var modStr in mods.Where(str => str.Length > 0))
                {
                    var modName = modStr.Split()[0];
                    var mod = Modification.Get(modName);
                    modComposition += mod.Composition;
                }

                if (ptmList.Count < 5)
                {
                    Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", scanNums[i], sequenceComp.Mass, modComposition.Mass, nMacthed[i], sequences[i], modsStr);
                }

                var compFromSeqAndMods = sequenceComp + modComposition;

                var expectedComposition = compositions[i];
                Assert.AreEqual(expectedComposition, compFromSeqAndMods.ToString(), "Composition Mismatch");

                ptmList.Add(new Tuple<int, double, double>(scanNums[i], sequenceComp.Mass, modComposition.Mass));
            }

            Console.WriteLine();
            Console.WriteLine("{0} results with PTMs / {1} total results passing threshold qValue < {2}", ptmList.Count, filterPassingResults, qValueThreshold);
            Assert.AreEqual(expectedPTMListCount, ptmList.Count, "Unexpected number of identifications with at least one PTM at qValue < {0}", qValueThreshold);
        }

        [Test]
        [TestCase(0.001, 12, 95)]
        [TestCase(0.01, 8, 98)]
        [TestCase(0.01, 12, 100)]
        [TestCase(0.01, 15, 100)]
        [TestCase(0.02, 12, 103)]
        public void TestParseMs1Ft(double qValueThreshold, int tolerancePpm, int expectedFeaturesWithId)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var resultFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt_IcTda.tsv");
            var resultFile = Utils.GetTestFile(methodName, resultFilePath);

            var ms1ftFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft");
            var ms1ftFile = Utils.GetTestFile(methodName, ms1ftFilePath);

            var resultParser = new MsPathFinderParser(resultFile.FullName);

            var idList = resultParser.GetIdList().TakeWhile(id => id.QValue <= qValueThreshold).OrderBy(id => id.Mass).ToList();
            var idMassList = idList.ConvertAll(id => id.Mass);
            var idFlag = new bool[idList.Count];

            var featureParser = new TsvFileParser(ms1ftFile.FullName);

            var minScan = featureParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            var maxScan = featureParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            var minCharge = featureParser.GetData("MinCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var maxCharge = featureParser.GetData("MaxCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var monoMass = featureParser.GetData("MonoMass").Select(Convert.ToDouble).ToArray();

            var numFeaturesWithId = 0;

            for (var i = 0; i < featureParser.NumData; i++)
            {
                var mass = monoMass[i];

                // Find Id
                var tolDa = new Tolerance(tolerancePpm).GetToleranceAsDa(mass, 1);
                var minMass = mass - tolDa;
                var maxMass = mass + tolDa;
                var index = idMassList.BinarySearch(mass);
                if (index < 0)
                {
                    index = ~index;
                }

                var matchedId = new List<MsPathFinderId>();

                // go down
                var curIndex = index - 1;
                while (curIndex >= 0)
                {
                    var curId = idList[curIndex];
                    if (curId.Mass < minMass)
                    {
                        break;
                    }

                    if (curId.Scan > minScan[i] && curId.Scan < maxScan[i]
                        && curId.Charge >= minCharge[i] && curId.Charge <= maxCharge[i])
                    {
                        matchedId.Add(curId);
                        idFlag[curIndex] = true;
                    }
                    --curIndex;
                }

                // go up
                curIndex = index;
                while (curIndex < idList.Count)
                {
                    var curId = idList[curIndex];
                    if (curId.Mass > maxMass)
                    {
                        break;
                    }

                    if (curId.Scan >= minScan[i] && curId.Scan <= maxScan[i]
                        && curId.Charge >= minCharge[i] && curId.Charge <= maxCharge[i])
                    {
                        matchedId.Add(curId);
                        idFlag[curIndex] = true;
                    }
                    ++curIndex;
                }

                if (matchedId.Count > 0)
                {
                    ++numFeaturesWithId;
                }
            }

            Console.WriteLine("Filtering on qValue < {0} and mass error < {1} ppm, find {2} features with an ID", qValueThreshold, tolerancePpm, numFeaturesWithId);
            Assert.AreEqual(expectedFeaturesWithId, numFeaturesWithId, "Unexpected number of features with an ID");
        }

        [Test]
        [Category("PNL_Domain")]
        public void TestFasta()
        {
            var db = new FastaDatabase(@"\\protoapps\UserData\Jungkap\Lewy\db\ID_005140_7A170668.fasta");
            Console.WriteLine(db.GetNumEntries());
        }

        [Test]
        [Category("PNL_Domain")]
        public void TestGeneratingMs1FeatureFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\QC_ShewIntact_1_19Jun15_Bane_14-09-01RZ.pbf";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var outFolderPath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, "Output");
            if (!Directory.Exists(outFolderPath))
            {
                Assert.Ignore("Skipping test {0} since folder not found: {1}", methodName, outFolderPath);
            }

            // string specFilePath = @"D:\MassSpecFiles\test\QC_Shew_Intact_4_01Jan15_Bane_C2-14-08-02RZ.raw";

            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 5000;
            const int maxThreads = 10;

            var param = new LcMsFeatureFinderInputParameters()
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
        public void TestFeatureMapGeneration()
        {
            Console.WriteLine("Testing Working");
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            var promexFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft");
            var promexFile = Utils.GetTestFile(methodName, promexFilePath);

            mFeatureMapPbfFile = pbfFile.FullName;
            mFeatureMapResultsFile = promexFile.FullName;

            var thread = new Thread(FeatureMapGeneration);

            thread.SetApartmentState(ApartmentState.STA);
            thread.Start();
            thread.Join();
        }

        private void FeatureMapGeneration()
        {
            var resultsFilePath = Path.Combine(Path.GetTempPath(), Path.GetFileNameWithoutExtension(mFeatureMapPbfFile) + "_FeatureMap.png");

            var map = new LcMsFeatureMap(PbfLcMsRun.GetLcMsRun(mFeatureMapPbfFile), mFeatureMapResultsFile, 2000, 50000);
            map.SaveImage(resultsFilePath);

            Console.WriteLine("Image saved to " + resultsFilePath);
        }

        [Test]
        [TestCase(6182.399, "508,510")]
        [TestCase(10868.732, "829")]
        [TestCase(17804.0293, "459,479,480,503,505")]
        public void TestProMexFilter(double massToFind, string expectedScanNumbers)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            var promexFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft");
            var promexFile = Utils.GetTestFile(methodName, promexFilePath);

            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName, 0, 0);

            var ms1Filter = new Ms1FtFilter(run, new Tolerance(10), promexFile.FullName);

            Console.WriteLine();

            var matchingScanNums = new SortedSet<int>();

            foreach (var item in ms1Filter.GetMatchingMs2ScanNums(massToFind))
            {
                matchingScanNums.Add(item);
            }

            var scanNumList = string.Join(",", matchingScanNums);

            Console.WriteLine("Scans with mass {0}:", massToFind);
            Console.WriteLine(scanNumList);

            var expectedScanNumList = expectedScanNumbers.Split(',');

            var matchCount = 0;
            foreach (var scanNumText in expectedScanNumList)
            {
                var scanNum = int.Parse(scanNumText);

                if (!matchingScanNums.Contains(scanNum))
                {
                    Assert.Fail("Did not find scan {0} for mass {1}", scanNum, massToFind);
                }

                matchCount++;
            }

            Assert.AreEqual(matchCount, matchingScanNums.Count, "Found extra matching scan numbers vs. what was expected");
        }

        [Test]
        [Category("PNL_Domain")]
        public void TestFeatureExampleForFigure()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFile = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf";
            // string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFile))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFile);
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var scorer = new LcMsFeatureLikelihood();
            var featureFinder = new LcMsPeakMatrix(run, scorer);
            var feature = featureFinder.GetLcMsPeakCluster(28061.6177, 20, 34, 7624, 7736);

            var resultsFilePath = Path.Combine(Path.GetTempPath(), Path.GetFileNameWithoutExtension(rawFile) + "_peaks.txt");
            var writer = new StreamWriter(resultsFilePath);

            writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n", "Scan", "Elution_Time", "Charge", "ID", "MZ", "Intensity", "Pearson_Correlation");

            var envelope = feature.TheoreticalEnvelope;
            foreach (var e in envelope.Isotopes)
            {
                Console.WriteLine(e.Ratio);
            }

            foreach (var env in feature.EnumerateEnvelopes())
            {
                var corr = env.PearsonCorrelation;
                for (var i = 0; i < envelope.Size; i++)
                {
                    var peak = env.Peaks[i];
                    if (peak == null)
                    {
                        continue;
                    }

                    writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n", env.ScanNum, run.GetElutionTime(env.ScanNum), env.Charge, i, peak.Mz, peak.Intensity, corr);
                }
            }
            writer.Close();

            Console.WriteLine("Results are in file " + resultsFilePath);
        }

        public static string[] TrainSetFileLists = {
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
