using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    class TestTopDownScoring
    {
        [Test]
        public void PrintAllScorers()
        {
            //            Console.WriteLine(Convert.ToDouble("0"));
            var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_SBEP.txt");
            scoringModel.PrintAllScores();
        }

        [Test]
        public void TestLikelihoodScorer()
        {
//            Console.WriteLine(Convert.ToDouble("0"));
            var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_Filtration_2.txt");
            Console.WriteLine("Score: {0}", scoringModel.GetScore(BaseIonType.Y, 0.99, 1200));
        }

        [Test]
        public void TestMatchedPeakCounter()
        {
            // Parameters
            var precursorIonTolerance = new Tolerance(15);
            var productIonTolerance = new Tolerance(15);

            var sw = new System.Diagnostics.Stopwatch();

            var aaSet = new AminoAcidSet();

            const string protAnnotation = "_.MFQQEVTITAPNGLHTRPAAQFVKEAKGFTSEITVTSNGKSASAKSLFKLQTLGLTQGTVVTISAEGEDEQKAVEHLVKLMAELE._";

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SBEP_STM_001_02272012_Aragon.raw";
            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

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
            var score = seqGraph.GetScore(charge, scorer);

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMostAbundantIsotopeMz(), ms2ScanNum, score);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestCorrMatchedPeakCounter()
        {
            // Parameters
            var precursorIonTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            var sw = new System.Diagnostics.Stopwatch();

            var aaSet = new AminoAcidSet();

            const string protAnnotation = "_.TMNITSKQMEITPAIRQHVADRLAKLEKWQTHLINPHIILSKEPQGFIADATINTPNGHLVASAKHEDMYTAINELINKLERQLNKVQHKGEAR._";

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SBEP_STM_001_02272012_Aragon.raw";
            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

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
            var score = seqGraph.GetScore(charge, scorer);

            Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMostAbundantIsotopeMz(), ms2ScanNum, score);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestMatchedPeakPostScorer()
        {
            // Parameters
            var productIonTolerance = new Tolerance(10);
            var scorer = new MatchedPeakPostScorer(productIonTolerance, 1, 10);
            var sw = new System.Diagnostics.Stopwatch();

            const int ms2ScanNum = 4658;
            var sequence = new Sequence("GYSIKDIIYQGEKSGVHNWQTLSGQNFYWHPDWLHIAEDLTGHKATASIQAEGTKATQNEAEQTIVKHLNKS", new AminoAcidSet());

            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            //const string specFilePath = @"D:\MassSpecFiles\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            Assert.True(spec != null);

            sw.Start();
            var score = scorer.ComputeScore(spec, sequence);

            Console.WriteLine("{0}\t{1}\t{2}", sequence, ms2ScanNum, score);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestJungkapScoring()
        {
            const string rawFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            var tolerance = new Tolerance(10);
            const int minCharge = 1;
            const int maxCharge = 15;

            var aminoAcidSet = new AminoAcidSet();
            var scorer = new MatchedPeakPostScorer(tolerance, minCharge, maxCharge);

            const string resultFileName = @"\\protoapps\UserData\Sangtae\TestData\IdFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTarget.tsv";
            var parser = new TsvFileParser(resultFileName);
            var scans = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var protSequences = parser.GetData("Sequence").ToArray();
            var modStrs = parser.GetData("Modifications").ToArray();
            var compositions = parser.GetData("Composition").Select(Composition.Parse).ToArray();

            const string outputFileName = @"\\protoapps\UserData\Sangtae\TestData\IdFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTarget_Rescored.tsv";
            using (var writer = new StreamWriter(outputFileName))
            {
                writer.WriteLine(string.Join(",", parser.GetHeaders()) + "\tScore");
                for (var i = 0; i < parser.NumData; i++)
                {
                    var scan = scans[i];
                    var protSequence = protSequences[i];
                    var modStr = modStrs[i];
                    //if (scan != 1765) continue;
                    var sequence = Sequence.CreateSequence(protSequence, modStr, aminoAcidSet);
                    //Console.WriteLine("{0}: {1} ? {2}", scan, sequence.Composition, compositions[i] - Composition.H2O);
                    Assert.True(sequence.Composition.Equals(compositions[i] - Composition.H2O));
                    var ms2Spec = run.GetSpectrum(scan) as ProductSpectrum;
                    Assert.True(ms2Spec != null);
                    var score = scorer.ComputeScore(ms2Spec, sequence);
                    writer.WriteLine("{0}\t{1}", parser.GetRows()[i], score);
                }
            }
            Console.WriteLine("Done");
        }
    }
}
