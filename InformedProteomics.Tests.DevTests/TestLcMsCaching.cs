using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests
{
    [TestFixture]
    internal class TestLcMsCaching
    {
        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        [TestCase(0.001, 108)]
        [TestCase(0.01, 113)]
        [TestCase(0.02, 116)]
        [TestCase(0.05, 130)]
        public void TestClusterCentricSearch(double qValueThreshold, int expectedNumCompositions)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var resultFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt_IcTda.tsv");
            var resultFile = Utils.GetTestFile(methodName, resultFilePath);

            var tsvReader = new TsvFileParser(resultFile.FullName);

            var ms2Scans = tsvReader.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var compositions = tsvReader.GetData("Composition").ToArray();
            var qValues = tsvReader.GetData("QValue").Select(Convert.ToDouble).ToArray();

            var compScanTable = new Dictionary<string, IList<int>>();
            for (var i = 0; i < qValues.Length; i++)
            {
                var qValue = qValues[i];
                if (qValue > qValueThreshold)
                {
                    break;
                }

                if (compScanTable.TryGetValue(compositions[i], out var scanNums))
                {
                    scanNums.Add(ms2Scans[i]);
                }
                else
                {
                    compScanTable.Add(compositions[i], new List<int> { ms2Scans[i] });
                }
            }

            Console.Write("NumCompositions: {0}", compScanTable.Keys.Count);

            Assert.AreEqual(expectedNumCompositions, compScanTable.Keys.Count);
        }

        [Test]
        [Category("PNL_Domain")]
        public void TestIsosFilter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var isosFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt_isos.csv");
            var isosfile = Utils.GetTestFile(methodName, isosFilePath);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName);
            var filter = new IsosFilter(run, new Tolerance(10), isosfile.FullName);

            const double massToFind = 944.08176;
            var matchingScanNums = filter.GetMatchingMs2ScanNums(massToFind).ToList();

            var scanNumList = string.Join(",", matchingScanNums);

            Console.WriteLine("Scans with mass {0}:", massToFind);
            Console.WriteLine(scanNumList);
        }

        [Test]
        [Category("Local_Testing")]
        public void FilteringEfficiencyQcShew()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, 1.4826, 1.4826);
            sw.Stop();

            Console.WriteLine("Reading run: {0:f4} sec", sw.Elapsed.TotalSeconds);

            const int minPrecursorCharge = 3;
            const int maxPrecursorCharge = 30;
            const int tolerancePpm = 10;
            var tolerance = new Tolerance(tolerancePpm);
            sw.Reset();
            sw.Start();
            var ms1BasedFilter = new Ms1IsotopeAndChargeCorrFilter(run, new Tolerance(10.0), minPrecursorCharge, maxPrecursorCharge, 3000, 50000, 0.7, 0.7, 0.7, 40);
            //var ms1BasedFilter = new Ms1IsotopeCorrFilter(run, minPrecursorCharge, maxPrecursorCharge, 15, 0.5, 40);

            sw.Stop();

            Console.WriteLine("Ms1 filter: {0:f4} sec", sw.Elapsed.TotalSeconds);

            ISequenceFilter ms1Filter = ms1BasedFilter;

            sw.Reset();
            sw.Start();
            const double minProteinMass = 3000.0;
            const double maxProteinMass = 30000.0;
            var minBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(minProteinMass);
            var maxBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(maxProteinMass);
            var numComparisons = 0L;
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                var mass = ProductScorerBasedOnDeconvolutedSpectra.GetMz(binNum);
                numComparisons += ms1Filter.GetMatchingMs2ScanNums(mass).Count();
            }
            sw.Stop();

            Console.WriteLine("Calculating #matches per bin: {0:f4} sec", sw.Elapsed.TotalSeconds);

            //const string prot =
            //    "ADVFHLGLTKAMLDGATLAIVPGDPERVKRIAELMDNATFLASHREYTSYLAYADGKPVVICSTGIGGPSTSIAVEELAQLGVNTFLRVGTTGAIQPHVNVGDVIVTQASVRLDGASLHFAPMEFPAVANFECTTAMVAACRDAGVEPHIGVTASSDTFYPGQERYDTVTGRVTRRFAGSMKEWQDMGVLNYEMESATLFTMCATQGWRAACVAGVIVNRTQQEIPDEATMKKTEVSAVSIVVAAAKKLLA";
            //var protMass = (new AminoAcidSet().GetComposition(prot) + Composition.H2O).Mass;
            //Console.WriteLine("************ScanNums: " + string.Join("\t", ms1Filter.GetMatchingMs2ScanNums(protMass)));

            const string resultFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\NoMod.tsv";
            if (!File.Exists(resultFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultFilePath);
            }

            var tsvReader = new TsvFileParser(resultFilePath);
            var scanNums = tsvReader.GetData("Scan(s)");
            var charges = tsvReader.GetData("Charge");
            var scores = tsvReader.GetData("E-value");
            var sequences = tsvReader.GetData("Peptide");

            //const string resultFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402_N30_C30.tsv";
            //var tsvReader = new TsvFileParser(resultFilePath);
            //var scanNums = tsvReader.GetData("ScanNum");
            //var charges = tsvReader.GetData("Charge");
            //var scores = tsvReader.GetData("Score");
            //var sequences = tsvReader.GetData("Sequence");

            var aaSet = new AminoAcidSet();

            var seqSet = new HashSet<string>();
            var allSeqSet = new HashSet<string>();
            var numUnfilteredSpecs = 0;
            var totalSpecs = 0;
            for (var i = 0; i < scores.Count; i++)
            {
                var score = Convert.ToDouble(scores[i]);
                if (score > 1E-4)
                {
                    continue;
                }
                //if (score < 10) continue;

                var scanNum = Convert.ToInt32(scanNums[i]);
                var charge = Convert.ToInt32(charges[i]);

                var sequence = SimpleStringProcessing.GetStringBetweenDots(sequences[i]);
                if (sequence == null || sequence.Contains("("))
                {
                    continue;
                }
                //var sequence = sequences[i];
                var composition = aaSet.GetComposition(sequence) + Composition.H2O;

                var precursorIon = new Ion(composition, charge);
                var isValid = run.GetSpectrum(scanNum) is ProductSpectrum spec && spec.IsolationWindow.Contains(precursorIon.GetMostAbundantIsotopeMz());
                if (!isValid)
                {
                    continue;
                }

                ++totalSpecs;

                var precursorScanNum = run.GetPrecursorScanNum(scanNum);
                var precursorSpec = run.GetSpectrum(precursorScanNum);
                var corr1 = precursorSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var nextScanNum = run.GetNextScanNum(scanNum, 1);
                var nextSpec = run.GetSpectrum(nextScanNum);
                var corr2 = nextSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var corr3 = ms1Filter.GetMatchingMs2ScanNums(composition.Mass).Contains(scanNum) ? 1 : 0;
                if (corr3 == 1)
                {
                    numUnfilteredSpecs++;
                    seqSet.Add(sequences[i]);
                }
                allSeqSet.Add(sequences[i]);

                var corrMax = new[] { corr1, corr2, corr3 }.Max();

                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", scanNum, precursorScanNum, corr1, nextScanNum, corr2, corr3, corrMax);
            }

            Console.WriteLine("TotalNumComparisons: {0}", numComparisons);
            Console.WriteLine("AverageNumComparisons: {0:f2}", numComparisons / (double)(maxBinNum - minBinNum + 1));
            Console.WriteLine("SuccessRate: {0:f2} {1} / {2}", numUnfilteredSpecs / (double)totalSpecs, numUnfilteredSpecs, totalSpecs);
            Console.WriteLine("NumUniqueSequences: {0:f2}, {1} / {2}", seqSet.Count / (double)allSeqSet.Count, seqSet.Count, allSeqSet.Count);

            Console.WriteLine("Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        [Category("Local_Testing")]
        public void FilteringEfficiency()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, 1.4826, 1.4826);
            sw.Stop();

            Console.WriteLine("Reading run: {0:f4} sec", sw.Elapsed.TotalSeconds);

            const int minPrecursorCharge = 3;
            const int maxPrecursorCharge = 30;
            const int tolerancePpm = 10;
            var tolerance = new Tolerance(tolerancePpm);
            sw.Reset();
            sw.Start();
            //var ms1BasedFilter = new Ms1BasedFilter(run, minPrecursorCharge, maxPrecursorCharge, tolerancePpm);
            //
            //var ms1BasedFilter = new Ms1IsotopeTopKFilter(run, minPrecursorCharge, maxPrecursorCharge, tolerancePpm, 20);
            //var ms1BasedFilter = new ProductScorerBasedOnDeconvolutedSpectra(run,
            //    minPrecursorCharge, maxPrecursorCharge,
            //    0, 0,
            //    600.0, 1800.0, new Tolerance(tolerancePpm), null);
            //ms1BasedFilter.CachePrecursorMatchesBinCentric();
            var ms1BasedFilter = new Ms1IsotopeAndChargeCorrFilter(run, new Tolerance(10.0), minPrecursorCharge, maxPrecursorCharge, 3000, 50000, 0.5, 0.5, 0.5, 40);
            //var ms1BasedFilter = new Ms1IsotopeCorrFilter(run, minPrecursorCharge, maxPrecursorCharge, 15, 0.5, 40);

            sw.Stop();

            Console.WriteLine("Ms1 filter: {0:f4} sec", sw.Elapsed.TotalSeconds);

            ISequenceFilter ms1Filter = ms1BasedFilter;

            sw.Reset();
            sw.Start();
            const double minProteinMass = 3000.0;
            const double maxProteinMass = 30000.0;
            var minBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(minProteinMass);
            var maxBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(maxProteinMass);
            var numComparisons = 0L;
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                var mass = ProductScorerBasedOnDeconvolutedSpectra.GetMz(binNum);
                numComparisons += ms1Filter.GetMatchingMs2ScanNums(mass).Count();
            }
            sw.Stop();

            Console.WriteLine("Calculating #matches per bin: {0:f4} sec", sw.Elapsed.TotalSeconds);

            const string resultFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon_4PTMs.icresult";
            if (!File.Exists(resultFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultFilePath);
            }

            var tsvReader = new TsvFileParser(resultFilePath);
            var compositions = tsvReader.GetData("Composition");
            var scanNums = tsvReader.GetData("ScanNum");
            var charges = tsvReader.GetData("Charge");
            var scores = tsvReader.GetData("Score");
            var qvalues = tsvReader.GetData("QValue");
            var sequences = tsvReader.GetData("Sequence");

            var sequenceCount = new Dictionary<string, int>();
            for (var i = 0; i < compositions.Count; i++)
            {
                if (qvalues != null)
                {
                    var qValue = Convert.ToDouble(qvalues[i]);
                    if (qValue > 0.01)
                    {
                        continue;
                    }
                }
                else
                {
                    var score = Convert.ToDouble(scores[i]);
                    if (score < 13)
                    {
                        continue;
                    }
                }
                var scanNum = Convert.ToInt32(scanNums[i]);
                var charge = Convert.ToInt32(charges[i]);
                var composition = Composition.Parse(compositions[i]);
                var precursorIon = new Ion(composition, charge);
                var isValid = run.GetSpectrum(scanNum) is ProductSpectrum spec && spec.IsolationWindow.Contains(precursorIon.GetMostAbundantIsotopeMz());
                if (!isValid)
                {
                    continue;
                }

                var sequence = sequences[i];
                if (sequenceCount.TryGetValue(sequence, out var count))
                {
                    sequenceCount[sequence] = count + 1;
                }
                else
                {
                    sequenceCount[sequence] = 1;
                }
            }
            //var sequences = tsvReader.GetData("Annotation");

            var seqSet = new HashSet<string>();
            var allSeqSet = new HashSet<string>();
            var numUnfilteredSpecs = 0;
            var totalSpecs = 0;
            for (var i = 0; i < compositions.Count; i++)
            {
                if (qvalues != null)
                {
                    var qValue = Convert.ToDouble(qvalues[i]);
                    if (qValue > 0.01)
                    {
                        continue;
                    }
                }
                else
                {
                    var score = Convert.ToDouble(scores[i]);
                    if (score < 13)
                    {
                        continue;
                    }
                }
                var scanNum = Convert.ToInt32(scanNums[i]);
                var charge = Convert.ToInt32(charges[i]);
                var composition = Composition.Parse(compositions[i]);
                var precursorIon = new Ion(composition, charge);
                var isValid = run.GetSpectrum(scanNum) is ProductSpectrum spec && spec.IsolationWindow.Contains(precursorIon.GetMostAbundantIsotopeMz());
                if (!isValid)
                {
                    continue;
                }

                ++totalSpecs;

                var precursorScanNum = run.GetPrecursorScanNum(scanNum);
                var precursorSpec = run.GetSpectrum(precursorScanNum);
                var corr1 = precursorSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var nextScanNum = run.GetNextScanNum(scanNum, 1);
                var nextSpec = run.GetSpectrum(nextScanNum);
                var corr2 = nextSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var corr3 = ms1Filter.GetMatchingMs2ScanNums(composition.Mass).Contains(scanNum) ? 1 : 0;
                if (corr3 == 1)
                {
                    numUnfilteredSpecs++;
                    seqSet.Add(sequences[i]);
                }
                allSeqSet.Add(sequences[i]);

                //var xic = run.GetFullPrecursorIonExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz(), tolerance);
                ////xic.Display();
                //var apexScanNum = xic.GetNearestApexScanNum(run.GetPrecursorScanNum(scanNum), false);
                //var apexSpec = run.GetSpectrum(apexScanNum);
                //var corr3 = apexSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var corrMax = new[] { corr1, corr2, corr3 }.Max();

                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}", scanNum, precursorScanNum, corr1, nextScanNum, corr2, corr3, corrMax, sequenceCount[sequences[i]]);
            }

            Console.WriteLine("TotalNumComparisons: {0}", numComparisons);
            Console.WriteLine("AverageNumComparisons: {0:f2}", numComparisons / (double)(maxBinNum - minBinNum + 1));
            Console.WriteLine("SuccessRate: {0:f2} {1} / {2}", numUnfilteredSpecs / (double)totalSpecs, numUnfilteredSpecs, totalSpecs);
            Console.WriteLine("NumUniqueSequences: {0:f2}, {1} / {2}", seqSet.Count / (double)allSeqSet.Count, seqSet.Count, allSeqSet.Count);

            Console.WriteLine("Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestFloatingPointRounding()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int numShifts = 36;
            const double value = 7655.9568537625;
            var converted = BitConverter.DoubleToInt64Bits(value);
            var rounded = (converted >> numShifts) << numShifts;
            var roundedDouble = BitConverter.Int64BitsToDouble(rounded);
            var roundedInt = (int)(rounded >> 32);
            Console.WriteLine("{0,25:E16}{1,23:X16}{2,23:X16}", value, converted, rounded);
            Console.WriteLine("{0}\t{1}", value, roundedDouble);
            Console.WriteLine("PPM error: {0}", (roundedDouble - value) / value * 1E6);
            Console.WriteLine("{0,16:X8}", roundedInt);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestPossibleSequenceMasses()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles\SBEP_STM_001_02272012_Aragon.raw";
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, 1.4826, 1.4826);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            //var ms1BasedFilter = new Ms1IsotopeCorrFilter(run, 3, 30, 15, 0.7, 1000);
            var ms1BasedFilter = new Ms1IsotopeAndChargeCorrFilter(run, new Tolerance(10));

            //var masses = ms1BasedFilter.GetPossibleSequenceMasses(1113);

            //var ms1BasedFilter = new Ms1IsotopeTopKFilter(run, 3, 30, 15);
            //var masses = ms1BasedFilter.GetPossibleSequenceMasses(2819, 20);
            //foreach (var m in masses)
            //{
            //    Console.WriteLine(m);
            //}
            sw.Stop();

            Console.WriteLine("Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestMs1Filtering()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string resultFilePath =
            //    @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrMatches_N30\SBEP_STM_001_02272012_Aragon.tsv";
                @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrMatches_N30\SBEP_STM_001_02272012_Aragon.decoy.icresult";
            if (!File.Exists(resultFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultFilePath);
            }

            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles\SBEP_STM_001_02272012_Aragon.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, 1.4826, 1.4826);

            //const int minPrecursorCharge = 3;
            //const int maxPrecursorCharge = 30;
            //const int tolerancePpm = 15;
            var tolerance = new Tolerance(15);

            //var ms1BasedFilter = new Ms1IsotopeCorrFilter(run, minPrecursorCharge, maxPrecursorCharge, tolerancePpm, 0.7, 40);
            ////var ms1BasedFilter = new Ms1IsotopeTopKFilter(run, minPrecursorCharge, maxPrecursorCharge, tolerancePpm, 20);
            //ISequenceFilter ms1Filter = ms1BasedFilter;

            var tsvReader = new TsvFileParser(resultFilePath);
            var compositions = tsvReader.GetData("Composition");
            var scanNums = tsvReader.GetData("ScanNum");
            var charges = tsvReader.GetData("Charge");
            var qValues = tsvReader.GetData("QValue");
            var scores = tsvReader.GetData("Score");

            //var sequences = tsvReader.GetData("Annotation");

            //var hist = new int[11];

            Console.WriteLine("ScanNum\tScore\tPrecursor\tNext\tSum\tNextIsotope\tLessCharge\tMoreCharge\tMax\tNumXicPeaks");
            for (var i = 0; i < compositions.Count; i++)
            {
                if (qValues != null)
                {
                    var qValue = Convert.ToDouble(qValues[i]);
                    if (qValue > 0.01)
                    {
                        continue;
                    }
                }

                var scanNum = Convert.ToInt32(scanNums[i]);
                var composition = Composition.Parse(compositions[i]);
                var charge = Convert.ToInt32(charges[i]);

                var precursorIon = new Ion(composition, charge);
                var isValid = run.GetSpectrum(scanNum) is ProductSpectrum spec && spec.IsolationWindow.Contains(precursorIon.GetMostAbundantIsotopeMz());
                if (!isValid)
                {
                    continue;
                }

                var score = Convert.ToDouble(scores[i]);

                var precursorScanNum = run.GetPrecursorScanNum(scanNum);
                var precursorSpec = run.GetSpectrum(precursorScanNum);
                var preIsotopeCorr = precursorSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var nextScanNum = run.GetNextScanNum(scanNum, 1);
                var nextSpec = run.GetSpectrum(nextScanNum);
                var nextIsotopeCorr = nextSpec.GetCorrScore(precursorIon, tolerance, 0.1);

                var xicMostAbundant = run.GetPrecursorExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz(), tolerance, scanNum);

                var apexScanNum = xicMostAbundant.GetApexScanNum();
                if (apexScanNum < run.MinLcScan)
                {
                    apexScanNum = scanNum;
                }
                //var sumSpec = run.GetSummedMs1Spectrum(apexScanNum);
                //                var apexIsotopeCorr = sumSpec.GetCorrScore(precursorIon, tolerance, 0.1);
                //                var corr3 = ms1Filter.GetMatchingMs2ScanNums(composition.Mass).Contains(scanNum) ? 1 : 0;

                var xicNextIsotope = run.GetPrecursorExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz() + Constants.C13MinusC12 / charge, tolerance, scanNum);

                var plusOneIsotopeCorr = xicMostAbundant.GetCorrelation(xicNextIsotope);

                var precursorIonChargeMinusOne = new Ion(composition, charge - 1);
                var xicChargeMinusOne = run.GetPrecursorExtractedIonChromatogram(precursorIonChargeMinusOne.GetMostAbundantIsotopeMz(), tolerance, scanNum);
                var chargeMinusOneCorr = xicMostAbundant.GetCorrelation(xicChargeMinusOne);

                var precursorIonChargePlusOne = new Ion(composition, charge + 1);
                var xicChargePlusOne = run.GetPrecursorExtractedIonChromatogram(precursorIonChargePlusOne.GetMostAbundantIsotopeMz(), tolerance, scanNum);
                var chargePlusOneCorr = xicMostAbundant.GetCorrelation(xicChargePlusOne);

                //var max = new[] {preIsotopeCorr, nextIsotopeCorr, apexIsotopeCorr, plusOneIsotopeCorr, chargeMinusOneCorr, chargePlusOneCorr}.Max();
                //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}",
                //    scanNum, score, preIsotopeCorr, nextIsotopeCorr, apexIsotopeCorr, plusOneIsotopeCorr, chargeMinusOneCorr, chargePlusOneCorr, max, xicMostAbundant.Count);
            }

            //Console.WriteLine("Histogram");
            //for (var i = 0; i < hist.Length; i++)
            //{
            //    Console.WriteLine("{0:f1}\t{1}", i / 10.0, hist[i]);
            //}
        }

        [Test]
        [Category("Local_Testing")]
        public void TestMs2Caching()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles\SBEP_STM_001_02272012_Aragon.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, 1.4826, 1.4826);

            //const int minPrecursorIonCharge = 3; // 3
            //const int maxPrecursorIonCharge = 30;// 67
            //const int minProductIonCharge = 1;
            //const int maxProductIonCharge = 10;

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            var runCache = new ProductScorerBasedOnDeconvolutedSpectra(run);
            runCache.DeconvoluteAllProductSpectra();
            sw.Stop();

            Console.WriteLine("Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestAveragine()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //for (var nominalMass = 1000; nominalMass <= 1000; nominalMass++)
            //{
            //    Console.WriteLine("{0}\t{1}", nominalMass,
            //        string.Join(",", Averagine.GetIsotopomerEnvelopeFromNominalMass(nominalMass).Envelope.Select(v => string.Format("{0:f3}", v))));
            //}
            for (var nominalMass = 1000; nominalMass <= 50000; nominalMass++)
            {
                var averagine = Averagine.GetIsotopomerEnvelopeFromNominalMass(nominalMass);
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void TestMs1Signature()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            const string resultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrMatches_N30";
            if (!File.Exists(resultPath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultPath);
            }

            foreach (var resultFilePath in Directory.GetFiles(resultPath, "*.tsv"))
            {
                Console.WriteLine(resultFilePath);
            }
        }

        [Test]
        [Category("PNL_Domain")]
        [Category("Long_Running")]
        public void TestNominalMassErrors()
        {

            const int MAX_RUNTIME_SECONDS = 60;

            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int minLength = 300;
            const int maxLength = 400;

            var sw = new System.Diagnostics.Stopwatch();

            var fastaFile = Utils.GetTestFile(methodName, Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"MSPathFinderT\ID_003962_71E1A1D4.fasta"));

            var db = new FastaDatabase(fastaFile.FullName);
            db.Read();
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            sw.Start();

            var hist = new long[11];
            var aaSet = new AminoAcidSet();
            foreach (var peptideAnnotationAndOffset in indexedDb.AnnotationsAndOffsetsNoEnzyme(minLength, maxLength))
            {
                ++numSequences;
                var annotation = peptideAnnotationAndOffset.Annotation;
                var sequenceStr = annotation.Substring(2, annotation.Length - 4);
                var sequenceComp = aaSet.GetComposition(sequenceStr);
                var mass = sequenceComp.Mass;
                var nominalMass = sequenceComp.NominalMass;
                var error = (int)Math.Round(mass * Constants.RescalingConstant) - nominalMass;
                var errorBin = error + hist.Length / 2;
                if (errorBin < 0)
                {
                    errorBin = 0;
                }

                if (errorBin >= hist.Length)
                {
                    errorBin = hist.Length - 1;
                }

                hist[errorBin]++;

                if (numSequences % 100 == 0 && sw.Elapsed.TotalSeconds > MAX_RUNTIME_SECONDS)
                {
                    break;
                }
            }

            Console.WriteLine("Sequence count: {0:N0}", numSequences);
            Console.WriteLine("{0,10}  {1,10}  {2,10}", "Bin ", "Count", "Fraction");
            for (var i = 0; i < hist.Length; i++)
            {
                Console.WriteLine("{0,10:F1}  {1,10:N0}  {2,10:F1}%", i - hist.Length / 2, hist[i], hist[i] / (double)numSequences * 100);
            }

            sw.Stop();

            Console.WriteLine("Elapsed Time: {0:F1} sec", sw.Elapsed.TotalSeconds);
        }
    }
}
