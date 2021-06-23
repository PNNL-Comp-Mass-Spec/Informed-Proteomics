using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    public class TestGeneratingFunction
    {
        // Ignore Spelling: Rescoring, const

        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        [TestCase(33, "MRHYEIVFMVHPDQSEQVPGMIERYTGAITEANGKIHRLEDWGRR")]
        [TestCase(446, "MKTFTATPETVTRDWFVVDADGKTLGRIATEIALRLRGKHKPEYTPHVDTGDYIIVINAEKVTVTGNKAQGKTYYSHSGFPGGIKQISFEKLQAHKPEMIIEKAVKGMLPKGPLGRAMFRKLKVYAGAEHNHAAQQPQVLDI")]
        public void TestGetScoreDistribution(int scanNum, string protSequence)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            if (!pbfFile.Exists)
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, pbfFile);
            }

            const string modStr = "";

            const int maxCharge = 20;
            const int minCharge = 1;
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);
            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName);

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
            //Console.WriteLine("{0}\t{1}", comparer.NumberOfBins, comparer.GetBinNumber(proteinMass));

            var stopwatch = Stopwatch.StartNew();
            var graphFactory = new ProteinScoringGraphFactory(comparer, aaSet);
            stopwatch.Stop();
            Console.WriteLine("edge generation elapsed time = {0:0.0} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

            var stopwatch2 = Stopwatch.StartNew();

            var sequence = Sequence.CreateSequence(protSequence, modStr, aaSet);
            var proteinMass = sequence.Mass + Composition.H2O.Mass;

            Console.WriteLine("Mass = {0}", proteinMass);

            var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;
            var deconvSpec = Deconvoluter.GetDeconvolutedSpectrum(spectrum, minCharge, maxCharge,
                isotopeOffsetTolerance, filteringWindowSize, tolerance, 0.7);

            stopwatch.Restart();

            var scorer = new CompositeScorerBasedOnDeconvolutedSpectrum(deconvSpec, spectrum, tolerance, comparer);
            var graph = graphFactory.CreateScoringGraph(scorer, proteinMass);
            stopwatch.Stop();
            Console.WriteLine("node generation elapsed time = {0:0.0} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

            stopwatch.Reset();
            stopwatch.Start();
            var gf = new GeneratingFunction(graph);
            gf.ComputeGeneratingFunction();
            //gf.ComputeGeneratingFunction(graph);
            stopwatch.Stop();
            Console.WriteLine("computing generation function = {0:0.0} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            var scoreDist = gf.GetScoreDistribution();

            Console.WriteLine("{0}-{1}", scoreDist.MinScore, scoreDist.MaxScore);

            Console.WriteLine("{0} : {1}", "score", "specEValue");

            for (var score = 15; score <= gf.MaximumScore; score++)
            {
                var specEvalue = gf.GetSpectralEValue(score);
                Console.WriteLine("{0} : {1}", score, specEvalue);
            }

            stopwatch2.Stop();
            Console.WriteLine("TOTAL computing generation function = {0:0.000} sec", stopwatch2.ElapsedMilliseconds / 1000.0d);
        }

        internal class TestMassBin : IMassBinning
        {
            internal TestMassBin()
            {
                MaxMass = 17000;
                MinMass = 0;
                NumberOfBins = GetBinNumber(MaxMass) - GetBinNumber(MinMass) + 1;
                Filtered = false;
            }

            public int GetBinNumber(double mass)
            {
                return Constants.GetBinNumHighPrecision(mass);
            }

            public double GetMass(int binNumber)
            {
                //throw new NotImplementedException();
                return binNumber / Constants.RescalingConstantHighPrecision;
            }

            public double GetMassStart(int binNumber)
            {
                return 0.5 * (GetMass(binNumber - 1) + GetMass(binNumber));
            }

            public double GetMassEnd(int binNumber)
            {
                return 0.5 * (GetMass(binNumber + 1) + GetMass(binNumber));
            }

            public double MaxMass { get; }
            public double MinMass { get; }
            public int NumberOfBins { get; }
            public bool Filtered { get; }
        }
        /*
        [Test]
        public void TestRescoring()
        {
            //const string specFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            const string specFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string sequence = "SGWYELSKSSNDQFKFVLKAGNGEVILTSELYTGKSGAMNGIESVQTNSPIEARYAKEVAKNDKPYFNLKAANHQIIGTSQMYSSTA";
            const int scanNum = 4084;
            const int charge = 7;
            var aminoAcidSet = new AminoAcidSet();
            var composition = aaSet.GetComposition(sequence) + Composition.H2O;

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            var informedScorer = new InformedTopDownScorer(run, aminoAcidSet, 1, 15, new Tolerance(10));
            var scores = informedScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, composition, charge, scanNum);
            Console.WriteLine(scores);
        }*/
    }
}
