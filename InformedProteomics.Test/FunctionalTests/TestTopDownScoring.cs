﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using DeconTools.Backend.ProcessingTasks.ChargeStateDeciders;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]

    class TestTopDownScoring
    {
        [Test]
        public void TestReadingMsDeconvFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";

            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const string filePath = @"H:\Research\QCShew_TopDown\Production\MsDeconvPlus\QC_Shew_Intact_26Sep14_Bane_C2Column3_msdeconv_plus.msalign";
            var parser = new MsDeconvFilter(run, new Tolerance(10), filePath);
        }
        /*
        [Test]
        public void PrintAllScorers()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //            Console.WriteLine(Convert.ToDouble("0"));
            const string filePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_SBEP.txt";            
            if (!File.Exists(filePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, filePath);
            }

            var scoringModel = new LikelihoodScoringModel(filePath);
            scoringModel.PrintAllScores();
        }

        [Test]
        public void TestLikelihoodScorer()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

//            Console.WriteLine(Convert.ToDouble("0"));

            const string filePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_Filtration_2.txt";
            if (!File.Exists(filePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, filePath);
            }

            var scoringModel = new LikelihoodScoringModel(filePath);
            Console.WriteLine("Score: {0}", scoringModel.GetScore(BaseIonType.Y, 0.99, 1200));
        }*/

        [Test]
        public void TestMatchedPeakCounter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Parameters
            var precursorIonTolerance = new Tolerance(15);
            var productIonTolerance = new Tolerance(15);

            var sw = new System.Diagnostics.Stopwatch();

            var aaSet = new AminoAcidSet();

            const string protAnnotation = "_.MFQQEVTITAPNGLHTRPAAQFVKEAKGFTSEITVTSNGKSASAKSLFKLQTLGLTQGTVVTISAEGEDEQKAVEHLVKLMAELE._";

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            Assert.NotNull(seqGraph, "Invalid sequence: {0}", protAnnotation);

            const string specFilePath = @"\\protoapps\UserData\Jungkap\Joshua\testData\SBEP_STM_001_02272012_Aragon.raw";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, 1.4826, 1.4826);

            sw.Start();
            var precursorFilter = new Ms1ContainsIonFilter(run, precursorIonTolerance);

            var seqCompositionArr = seqGraph.GetSequenceCompositions();
            Console.WriteLine("Length: {0}\tNumCompositions: {1}", protAnnotation.Length - 4, seqCompositionArr.Length);

            const int charge = 6;
            const int modIndex = 0;
            const int ms2ScanNum = 4448;

            var seqComposition = seqCompositionArr[modIndex];
            var peptideComposition = seqComposition + Composition.H2O;
            peptideComposition.GetIsotopomerEnvelopeRelativeIntensities();

            Console.WriteLine("Composition: {0}, AveragineMass: {1}", seqComposition, seqComposition.Mass);
            seqGraph.SetSink(modIndex);

            var precursorIon = new Ion(peptideComposition, charge);

            Assert.True(precursorFilter.IsValid(precursorIon, ms2ScanNum));

            var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            Assert.True(spec != null);

            var scorer = new MatchedPeakCounter(spec, productIonTolerance, 1, 10);
            var score = seqGraph.GetFragmentScore(scorer);

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMostAbundantIsotopeMz(), ms2ScanNum, score);

            sw.Stop();

            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestCorrMatchedPeakCounter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Parameters
            var precursorIonTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            var sw = new System.Diagnostics.Stopwatch();

            var aaSet = new AminoAcidSet();

            const string protAnnotation = "_.TMNITSKQMEITPAIRQHVADRLAKLEKWQTHLINPHIILSKEPQGFIADATINTPNGHLVASAKHEDMYTAINELINKLERQLNKVQHKGEAR._";

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            Assert.NotNull(seqGraph, "Invalid sequence: {0}", protAnnotation);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SBEP_STM_001_02272012_Aragon.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, 1.4826, 1.4826);

            sw.Start();
            var precursorFilter = new Ms1ContainsIonFilter(run, precursorIonTolerance);

            var seqCompositionArr = seqGraph.GetSequenceCompositions();
            Console.WriteLine("Length: {0}\tNumCompositions: {1}", protAnnotation.Length - 4, seqCompositionArr.Length);

            const int charge = 9;
            const int modIndex = 0;
            const int ms2ScanNum = 3633;

            var seqComposition = seqCompositionArr[modIndex];
            var peptideComposition = seqComposition + Composition.H2O;
            peptideComposition.GetIsotopomerEnvelopeRelativeIntensities();

            Console.WriteLine("Composition: {0}, AveragineMass: {1}", seqComposition, seqComposition.Mass);
            seqGraph.SetSink(modIndex);

            var precursorIon = new Ion(peptideComposition, charge);

            Assert.True(precursorFilter.IsValid(precursorIon, ms2ScanNum));

            var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            Assert.True(spec != null);

            //var scorer = new MatchedPeakCounter(spec, productIonTolerance, 1, 10);
            var scorer = new CorrMatchedPeakCounter(spec, productIonTolerance, 1, 10);
            var score = seqGraph.GetFragmentScore(scorer);

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMostAbundantIsotopeMz(), ms2ScanNum, score);

            sw.Stop();

            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestMatchedPeakPostScorer()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Parameters
            var productIonTolerance = new Tolerance(10);
            var scorer = new MatchedPeakPostScorer(productIonTolerance, 1, 10);
            var sw = new System.Diagnostics.Stopwatch();

            const int ms2ScanNum = 4658;
            var sequence = new Sequence("GYSIKDIIYQGEKSGVHNWQTLSGQNFYWHPDWLHIAEDLTGHKATASIQAEGTKATQNEAEQTIVKHLNKS", new AminoAcidSet());

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            //const string specFilePath = @"D:\MassSpecFiles\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            Assert.True(spec != null);

            sw.Start();
            var score = scorer.ComputeScore(spec, sequence);

            Console.WriteLine("{0}\t{1}\t{2}", sequence, ms2ScanNum, score);

            sw.Stop();

            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestCompositeScoring()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string rawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            const string rawFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            // Configure amino acid set
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                oxM,
                acetylN
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            var comparer = new FilteredProteinMassBinning(aaSet, 50000, 28);

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);
            const int minCharge = 1;
            const int maxCharge = 20;
            var graphFactory = new ProteinScoringGraphFactory(comparer, aaSet);
            var aminoAcidSet = new AminoAcidSet();
            //var scorer = new MatchedPeakPostScorer(tolerance, minCharge, maxCharge);
            var scorer = new InformedTopDownScorer(run, aminoAcidSet, minCharge, maxCharge, tolerance);

            var fileExt = new string[] {"IcTarget", "IcDecoy"};
            foreach (var ext in fileExt)
            {
                var resultFileName = string.Format(@"D:\MassSpecFiles\training\Rescoring\QC_Shew_Intact_26Sep14_Bane_C2Column3_{0}.tsv", ext);
                var parser = new TsvFileParser(resultFileName);
                var scans = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
                var charges = parser.GetData("Charge").Select(s => Convert.ToInt32(s)).ToArray();
                var protSequences = parser.GetData("Sequence").ToArray();
                var modStrs = parser.GetData("Modifications").ToArray();
                var compositions = parser.GetData("Composition").Select(Composition.Parse).ToArray();
                var protMass = parser.GetData("Mass").Select(s => Convert.ToDouble(s)).ToArray();
                var outputFileName = string.Format(@"D:\MassSpecFiles\training\Rescoring\QC_Shew_Intact_26Sep14_Bane_C2Column3_{0}_Rescored.tsv", ext);

                using (var writer = new StreamWriter(outputFileName))
                {
                    writer.WriteLine(string.Join("\t", parser.GetHeaders().ToArray(), 0, 15) + "\tScore\tEValue");

                    var lines = new string[parser.NumData];
                    
                    //for (var i = 0; i < parser.NumData; i++)
                    Parallel.For(0, parser.NumData, i =>
                    {
                        var scan = scans[i];
                        var charge = charges[i];
                        var protSequence = protSequences[i];
                        var modStr = modStrs[i];
                        var sequence = Sequence.CreateSequence(protSequence, modStr, aminoAcidSet);
                        Assert.True(sequence.Composition.Equals(compositions[i] - Composition.H2O));
                        var ms2Spec = run.GetSpectrum(scan) as ProductSpectrum;
                        Assert.True(ms2Spec != null);
                        var scores = scorer.GetScores(sequence, charge, scan);

                        var deconvSpec = Deconvoluter.GetDeconvolutedSpectrum(ms2Spec, minCharge, maxCharge,
                            isotopeOffsetTolerance, filteringWindowSize, tolerance, 0.7);

                        var deconvScorer = new CompositeScorerBasedOnDeconvolutedSpectrum(deconvSpec, ms2Spec, tolerance,
                            comparer);
                        var graph = graphFactory.CreateScoringGraph(deconvScorer, protMass[i]);

                        var gf = new GeneratingFunction(graph);
                        gf.ComputeGeneratingFunction();

                        var specEvalue = gf.GetSpectralEValue(scores.Score);

                        var rowStr = parser.GetRows()[i];
                        var items = rowStr.Split('\t').ToArray();
                        var newRowStr = string.Join("\t", items, 0, 15);

                        //writer.WriteLine("{0}\t{1}\t{2}", newRowStr, scores.Score, specEvalue);
                        lock (lines)
                        {
                            lines[i] = string.Format("{0}\t{1}\t{2}", newRowStr, scores.Score, specEvalue);    
                        }
                        //Console.WriteLine("{0}\t{1}\t{2}", items[0], scores.Score, specEvalue);
                    });

                    foreach (var line in lines) writer.WriteLine(line);
                }
                Console.WriteLine("Done");
            }
        }

        [Test]
        public void RecomputeFdr()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string targetResultPath = @"D:\MassSpecFiles\training\Rescoring\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTarget_Rescored.tsv";
            const string decoyResultPath = @"D:\MassSpecFiles\training\Rescoring\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcDecoy_Rescored.tsv";
            const string tdaResultPath = @"D:\MassSpecFiles\training\Rescoring\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda_Rescored.tsv";

            //const string targetResultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.icresult";
            //const string decoyResultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.decoy.icresult";
            var fdrCalculator = new FdrCalculator(targetResultPath, decoyResultPath);
            if (fdrCalculator.HasError())
            {
                throw new Exception(@"Error computing FDR: " + fdrCalculator.ErrorMessage);
            }

            fdrCalculator.WriteTo(tdaResultPath);
            Console.WriteLine("Done");
        }

    }
}
