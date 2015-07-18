using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestYufengData
    {
        public const string TestRawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\yufeng_column_test2.raw";

        [Test]
        public void AddMostAbundantIsotopePeakIntensity()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\QC_ShewIntact_40K_LongSeparation_1_141016155143.raw";

            if (!File.Exists(rawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + rawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const string resultFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\QC_ShewIntact_40K_LongSeparation_1_141016155143_IcTda.tsv";

            var parser = new TsvFileParser(resultFilePath);
            var compositions = parser.GetData("Composition").Select(Composition.Parse).ToArray();
            var scanNums = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var charges = parser.GetData("Charge").Select(s => Convert.ToInt32(s)).ToArray();
            var precursorIntensities = new double[parser.NumData];
            var tolerance = new Tolerance(10);
            for (var i = 0; i < parser.NumData; i++)
            {
                var scanNum = scanNums[i];
                var composition = compositions[i];
                var charge = charges[i];
                var precursorIon = new Ion(composition, charge);

                var precursorScanNum = run.GetPrecursorScanNum(scanNum);
                var precursorSpec = run.GetSpectrum(precursorScanNum);
                var isotopePeaks = precursorSpec.GetAllIsotopePeaks(precursorIon, tolerance, 0.1);
                if (isotopePeaks != null)
                {
                    var maxIntensity = 0.0;
                    for (var j = 0; j < isotopePeaks.Length; j++)
                    {
                        if (isotopePeaks[j] != null && isotopePeaks[j].Intensity > maxIntensity)
                            maxIntensity = isotopePeaks[j].Intensity;
                    }
                    precursorIntensities[i] = maxIntensity;
                }
            }

            // Writing
            const string newResultFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\QC_ShewIntact_40K_LongSeparation_1_141016155143_IcTdaWithIntensities.tsv";
            using (var writer = new StreamWriter(newResultFilePath))
            {
                writer.WriteLine(string.Join("\t", parser.GetHeaders())+"\t"+"PrecursorIntensity");
                for (var i = 0; i < parser.NumData; i++)
                {
                    writer.WriteLine(parser.GetRows()[i]+"\t"+precursorIntensities[i]);
                }
            }
            Console.WriteLine("Done");
        }

        [Test]
        public void Test43KProtein()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Configure amino acid set
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var dethiomethylM = new SearchModification(Modification.Dethiomethyl, 'M', SequenceLocation.Everywhere, false);
            var deamidatedN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            var deamidatedQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);
            var pyroCarbamidomethylC = new SearchModification(Modification.PyroCarbamidomethyl, 'C',
                SequenceLocation.ProteinNTerm, false);
            var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                oxM,
                dethiomethylM,
                acetylN,
                //phosphoS,
                //phosphoT,
                //phosphoY,
                deamidatedN,
//                deamidatedQ,
                glutathioneC,
                pyroCarbamidomethylC,
                nitrosylC,
                nethylmaleimideC
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
//            var aaSet = new AminoAcidSet();

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath);
            const string protSequence =
                "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(aaSet, AminoAcid.ProteinNTerm, protSequence, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return;

            var ms1Filter = new SimpleMs1Filter();
            var ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(run);
            foreach(var ms2ScanNum in Ms2ScanNums) ms2ScorerFactory.GetScorer(ms2ScanNum);

            for (var numNTermCleavages = 0; numNTermCleavages <= 0; numNTermCleavages++)
            {
                if (numNTermCleavages > 0) seqGraph.CleaveNTerm();
                var numProteoforms = seqGraph.GetNumProteoforms();
                var modCombs = seqGraph.GetModificationCombinations();
                for (var modIndex = 0; modIndex < numProteoforms; modIndex++)
                {
                    seqGraph.SetSink(modIndex);
                    var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    var sequenceMass = protCompositionWithH2O.Mass;
                    var modCombinations = modCombs[modIndex];

                    foreach (var ms2ScanNum in ms1Filter.GetMatchingMs2ScanNums(sequenceMass))
                    {
                        var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                        if (spec == null) continue;
                        var charge =
                            (int)
                                Math.Round(sequenceMass /
                                           (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));
                        var scorer = ms2ScorerFactory.GetMs2Scorer(ms2ScanNum);
                        var score = seqGraph.GetFragmentScore(scorer);
                        if (score <= 3) continue;

                        var precursorIon = new Ion(protCompositionWithH2O, charge);
                        var sequence = protSequence.Substring(numNTermCleavages);
                        var pre = numNTermCleavages == 0 ? annotation[0] : annotation[numNTermCleavages + 1];
                        var post = annotation[annotation.Length - 1];

                        Console.WriteLine("{0}.{1}.{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}", pre, sequence, post, ms2ScanNum, modCombinations, 
                            precursorIon.GetMostAbundantIsotopeMz(), precursorIon.Charge, precursorIon.Composition.Mass, score);
                    }
                }
            }
        }

        private static readonly IEnumerable<int> Ms2ScanNums = new[] { 46454, 46475, 46484, 46506, 46562, 46661 };
        private class SimpleMs1Filter : ISequenceFilter
        {
            public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
            {
                return Ms2ScanNums;
            }
        }

        [Test]
        public void TestDeconvolution()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            const int minScanNum = 46454;   // 635.43
            const int maxScanNum = 46661;   // 638.90
            const int MAX_POINTS = 50;

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath) as PbfLcMsRun;
            if (run == null) return;
            var summedSpec = run.GetSummedMs1Spectrum(minScanNum, maxScanNum);
            summedSpec.FilterNoise(50.0);
            // summedSpec.Display(MAX_POINTS);

            var deconvoluted = ProductScorerBasedOnDeconvolutedSpectra.GetDeconvolutedSpectrum(summedSpec, 2, 45, new Tolerance(10), 0.9, 2);
            deconvoluted.Display(MAX_POINTS);
        }

        [Test]
        public void TestSumMs1Spectra()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            const int minScanNum = 46454;
            const int maxScanNum = 46661;
            const int MAX_POINTS = 50;

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath) as PbfLcMsRun;
            if (run == null) return;
            var summedSpec = run.GetSummedMs1Spectrum(minScanNum, maxScanNum);

            summedSpec.Display(MAX_POINTS);
        }

        [Test]
        public void TestRunningTimeSummingSpectra()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath, 1.4826, 1.4826) as PbfLcMsRun;

            var sw = new Stopwatch();
            sw.Start();

            const int windowSize = 5;
            foreach (var scanNum in run.GetScanNumbers(1))
            {
                //var spec = run.GetSpectrum(scanNum);
                var spec = run.GetSummedMs1Spectrum(Math.Max(scanNum - windowSize, run.MinLcScan),
                    Math.Min(scanNum + windowSize, run.MaxLcScan));
            }
            sw.Stop();
            
            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestSumIsoProfilesAcrossDifferentCharges()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath) as PbfLcMsRun;

            //var spec = run.GetSpectrum(46452); // 635.37
            var spec = run.GetSummedMs1Spectrum(46437, 46466);
            var tolerance = new Tolerance(10);

            const string protSequence =
                "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            //const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(new AminoAcidSet(), AminoAcid.ProteinNTerm, protSequence, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return;
            seqGraph.SetSink(0);
            var neutral = seqGraph.GetSinkSequenceCompositionWithH2O();

            var theoProfile = neutral.GetIsotopomerEnvelopeRelativeIntensities();
            var expProfile = new double[theoProfile.Length];
            for (var charge = 22; charge <= 45; charge++)
            {
                var ion = new Ion(neutral, charge);
                var isotopePeaks = spec.GetAllIsotopePeaks(ion, tolerance, 0.1);
                if (isotopePeaks == null) continue;
                Assert.True(isotopePeaks.Length == theoProfile.Length);
                for (var i = 0; i < isotopePeaks.Length; i++)
                {
                    if (isotopePeaks[i] != null) expProfile[i] += isotopePeaks[i].Intensity;
                }
            }
            for (var i = 0; i < theoProfile.Length; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}", neutral.GetIsotopeMass(i), theoProfile[i], expProfile[i]/expProfile.Max());
            }
            Console.WriteLine("Corr: " + FitScoreCalculator.GetPearsonCorrelation(theoProfile, expProfile));
        }

        [Test]
        public void TestSmartIsoWindowSumming()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            const string protSequence =
                "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            //const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(new AminoAcidSet(), AminoAcid.ProteinNTerm, protSequence, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return;
            seqGraph.SetSink(0);
            var neutral = seqGraph.GetSinkSequenceCompositionWithH2O();
            var ion = new Ion(neutral, 43);
            var tolerance = new Tolerance(10);

            const int targetMs2ScanNum = 46562;
            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath) as PbfLcMsRun;
            var ms2Spec = run.GetSpectrum(targetMs2ScanNum) as ProductSpectrum;
            Assert.True(ms2Spec != null);
            var isoWindow = ms2Spec.IsolationWindow;
            //var prevScanNum = run.GetPrevScanNum(targetMs2ScanNum, 1);
            //var nextScanNum = run.GetNextScanNum(targetMs2ScanNum, 1);
            var summedSpec = run.GetSummedMs1Spectrum(targetMs2ScanNum, 2.5);
            //var windowSpec = summedSpec.GetPeakListWithin(isoWindow.MinMz, isoWindow.MaxMz);
            Console.WriteLine("Corr: " + summedSpec.GetCorrScore(ion, tolerance));
        }

        public void TestCorrelation()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            Console.WriteLine(@"Test not implemented: " + methodName);
        }

        [Test]
        public void GetIsoProfile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string protSequence =
                "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            //const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(new AminoAcidSet(), AminoAcid.ProteinNTerm, protSequence, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return;
            seqGraph.SetSink(0);
            var neutral = seqGraph.GetSinkSequenceCompositionWithH2O() - Composition.Hydrogen;
            //Console.WriteLine(neutral);

            for (var charge = 22; charge <= 60; charge++)
            {
                var ion = new Ion(neutral, charge);
                Console.WriteLine("{0}\t{1}", charge, ion.GetMostAbundantIsotopeMz());
            }

            var ion27 = new Ion(neutral, 29);
            var isotopes = ion27.GetIsotopes(0.1); 
            foreach (var isotope in isotopes)
            {
                Console.WriteLine("{0}\t{1}", ion27.GetIsotopeMz(isotope.Index), isotope.Ratio);
            }                
        }

        [Test]
        public void TestGetNumBins()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var comparer = new MzComparerWithBinning(26);
            const double minMz = 5000.0; // 600.0
            const double maxMz = 10000.0;    // 2000.0
            var minBinNum = comparer.GetBinNumber(minMz);
            var maxBinNum = comparer.GetBinNumber(maxMz);
            var numBins = maxBinNum - minBinNum;
            Console.WriteLine(@"NumBins: " + numBins);

            Assert.IsTrue(numBins == 16384);
        }        

        [Test]
        public void TestGeneratingXicsOfAllCharges()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }
            
            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath, 0.0, 0.0);
            var comparer = new MzComparerWithBinning(27);
            const string protSequence =
                "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            //const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(new AminoAcidSet(), AminoAcid.ProteinNTerm, protSequence, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return;
            seqGraph.SetSink(0);
            var neutral = seqGraph.GetSinkSequenceCompositionWithH2O() - Composition.Hydrogen;
            var proteinMass = neutral.Mass;
            var isoEnv = Averagine.GetIsotopomerEnvelope(proteinMass);

            const bool SHOW_ALL_SCANS = false;
            var targetColIndex = 0;

            if (SHOW_ALL_SCANS)           
                Console.WriteLine("Charge\t" + string.Join("\t", run.GetScanNumbers(1)));
            else
            {
                // Just display data for scan 161
                Console.WriteLine("Charge\t161");
                foreach (var scanNumber in run.GetScanNumbers(1))
                {
                    if (scanNumber == 161)
                        break;
                    targetColIndex++;
                }

            }

            const int minCharge = 2;
            const int maxCharge = 60;
            for (var charge = minCharge; charge <= maxCharge; charge++)
            {
                var ion = new Ion(neutral, charge);
                var mostAbundantIsotopeMz = ion.GetIsotopeMz(isoEnv.MostAbundantIsotopeIndex);
                //var secondMostAbundantIsotopeMz = ion.GetIsotopeMz(isoEnv.MostAbundantIsotopeIndex + 1);
                var binNum = comparer.GetBinNumber(mostAbundantIsotopeMz);
                var mzStart = comparer.GetMzStart(binNum);
                var mzEnd = comparer.GetMzEnd(binNum);

                var xic = run.GetFullPrecursorIonExtractedIonChromatogram(mzStart, mzEnd);
                Console.Write(charge+"\t");

                
                if (SHOW_ALL_SCANS)
                    Console.WriteLine(string.Join("\t", xic.Select(p => p.Intensity)));
                else
                    Console.WriteLine(xic[targetColIndex].Intensity);
            }
        }
      
        [Test]
        public void TestAbpSumMs1Spectra()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";

            if (!File.Exists(specFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + specFilePath);
                return;
            }

            const int minScanNum = 5657;
            const int maxScanNum = 5699;
            const int MAX_POINTS = 50;

            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
            if (run == null) return;
            var summedSpec = run.GetSummedMs1Spectrum(minScanNum, maxScanNum);
            var peakList = summedSpec.GetPeakListWithin(1180.0, 1192.0);
            var filteredPeakList = new List<Peak>();
            PeakListUtils.FilterNoise(peakList, ref filteredPeakList);
            new Spectrum(filteredPeakList, 0).Display(MAX_POINTS);
        }

        [Test]
        public void TestIsoProfile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string sequence = "MWYMISAQDVENSLEKRLAARPAHLARLQELADEGRLLVAGPHPAIDSENPGDAGFSGSLVVADFDSLATAQAWADADPYFAAGVYQSVVVKPFKRVLP";
            var aaSet = new AminoAcidSet();
            var comp = aaSet.GetComposition(sequence) + Composition.H2O;
            var ion = new Ion(comp, 9);
            foreach (var i in ion.GetIsotopes(0.1))
            {
                Console.WriteLine(ion.GetIsotopeMz(i.Index)+"\t"+i.Ratio);
            }
        }

        [Test]
        public void TestSumMs2Spectra()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TestYufengData\NewQC_LongSep_29Sep14_141001104925.raw";

            if (!File.Exists(specFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + specFilePath);
                return;
            }

            const int minScanNum = 1289;
            //const int maxScanNum = 1389;
            const int minCharge = 6;
            //const int maxCharge = 6;
            const string sequence = "EIRGYRPPEPYKGKGVRYDDEEVRRKEAKKK";
            var aaSet = new AminoAcidSet();
            
            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
            
            var scorer = new InformedTopDownScorer(run, aaSet, 1, minCharge - 1, new Tolerance(10));
            scorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm,
                Composition.Parse("C(166) H(270) N(52) O(49) S(0)"), minCharge, minScanNum);
        }
    }
}
