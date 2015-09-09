using System;
using System.Collections;
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
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using NUnit.Framework;
using Peak = InformedProteomics.Backend.Data.Spectrometry.Peak;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestGeneratingFunction
    {
        
        [Test]
        public void TestGetScoreDistribution()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            //const string rawFile = @"D:\MassSpecFiles\training\raw\yufeng_column_test2.pbf";

            const int scanNum = 5109;
            const string protSequence = "ENIAVVDMGAVFEQLPQREQIMQSLKSEFGDRMSEVQKMQEEMRSLMEKQQRDGALMNDTQKTELVRKMEALKSEYQLKGKALDEDLRRRQGEEQNKLLVKVQKAINTIAEKEKYDLVLQRGAVIYVKPNADISGKVVEALSKGK";
            
            //const string protSequence =
            //    "AIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVGLHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTVTSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVGIGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGSAAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEANQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDLKSLKELLKDQEGAVALKIVRGKSMLYLVLR";
            
            const string modStr = "";
            
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
            var sequence = Sequence.CreateSequence(protSequence, modStr, aaSet);
            var proteinMass = sequence.Mass + Composition.H2O.Mass;

            Console.WriteLine("Mass = {0}", proteinMass);
            

            const int maxCharge = 20;
            const int minCharge = 1;
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

            if (!File.Exists(rawFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFile);
            }

            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spectrum, minCharge, maxCharge,
                isotopeOffsetTolerance, filteringWindowSize, tolerance, 0.7);

            var peakList = deconvolutedPeaks.Select(dp => new Peak(dp.Mass, dp.Intensity)).ToArray();
            var productSpec = new ProductSpectrum(peakList, spectrum.ScanNum)
            {
                MsLevel = spectrum.MsLevel,
                ActivationMethod = spectrum.ActivationMethod,
                IsolationWindow = spectrum.IsolationWindow
            };

            //var comparer = new FilteredProteinMassBinning(aaSet, 50000, 28);
            //var comparer = new ProteinMassBinning(50, 50000);
            var comparer = new TestMassBin();
            Console.WriteLine("{0}\t{1}",comparer.NumberOfBins, comparer.GetBinNumber(proteinMass));
            var graphFactory = new ProteinScoringGraphFactory(comparer, aaSet);

            Console.WriteLine(proteinMass);
            var stopwatch = Stopwatch.StartNew();
            var graph = graphFactory.CreateScoringGraph(productSpec, proteinMass);
            stopwatch.Stop();
            Console.WriteLine(@"graph generation elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            
            stopwatch.Reset();
            stopwatch.Start();
            var gf = new GeneratingFunction(comparer.NumberOfBins);
            gf.ComputeGeneratingFunction(graph);
            stopwatch.Stop();
            Console.WriteLine(@"computing generation function = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            for (var score = 1; score <= gf.MaximumScore; score++)
            {
                var specEvalue = gf.GetSpectralEValue(score);
                Console.WriteLine("{0} : {1}", score, specEvalue);                                
            }


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
                return binNumber/Constants.RescalingConstantHighPrecision;
            }

            public double GetMassStart(int binNumber)
            {
                return 0.5*(GetMass(binNumber - 1) + GetMass(binNumber));

            }

            public double GetMassEnd(int binNumber)
            {
                return 0.5 * (GetMass(binNumber + 1) + GetMass(binNumber));
            }

            public double MaxMass { get; private set; }
            public double MinMass { get; private set; }
            public int NumberOfBins { get; private set; }
            public bool Filtered { get; private set; }
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
            var aaSet = new AminoAcidSet();
            var composition = aaSet.GetComposition(sequence) + Composition.H2O;

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            var informedScorer = new InformedTopDownScorer(run, aaSet, 1, 15, new Tolerance(10));
            var scores = informedScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, composition, charge, scanNum);
            Console.WriteLine(scores);
        }*/
    }



}
