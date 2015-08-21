using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;
using PNNLOmics.Data;
using Peak = InformedProteomics.Backend.Data.Spectrometry.Peak;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestGeneratingFunction
    {
        [Test]
        public void Test()
        {
            var massList = AminoAcid.StandardAminoAcidArr.Select(aa => aa.Mass).ToList();
            var MaxMass = 50000;
            var fineNodes = new BitArray(Constants.GetBinNumHighPrecision(MaxMass));
            fineNodes[0] = true;
            
            var effectiveBinCounter = 0;
            for (var fineBinIdx = 0; fineBinIdx < fineNodes.Length; fineBinIdx++)
            {
                if (!fineNodes[fineBinIdx]) continue;

                var fineNodeMass = fineBinIdx/Constants.RescalingConstantHighPrecision;
                foreach (var m in massList)
                {
                    var validFineNodeIndex = Constants.GetBinNumHighPrecision(fineNodeMass + m);
                    if (validFineNodeIndex >= fineNodes.Length) break;
                    fineNodes[validFineNodeIndex] = true;
                }
            }

            Console.WriteLine(fineNodes.Length);
            var n = 0;
            for(var i = 0; i < fineNodes.Length; i++) if (fineNodes[i]) n++;
            Console.WriteLine(n);

        }
        
        [Test]
        public void TestGetScoreDistribution()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            
            const int scanNum = 6745;
            const string protSequence = "EKQRPDFIYPQGSMLMHSPLVLNKADMYGFFLKGKLANLQRSIDNTLNQVAGSAMYFKVLSPYVLTTFTQIEKAYSDYPTDRAKGWIQETDIITWVMVGRQDNSSSTKISHVYFQPLHIWVNDAMALINGRELFGYPKYLCEYTMPAAGQPLTRLSLAAKSFLHFSPDTELAMHPLLEVNCNAKPEAQLSTVEAITQTWKLFKEQTDFIPELDKLGEEQLFHLLFKPAIDQVFLKQFPDASGQKAVYQAITASPAKVNKVHSVALLENDFVATLFDNASFPLKDTLGVELGEQHVLLPYHVNFDFEVPPGEVLVDNSEIKKQKIAILGGGVAAMTAACYLTDKPGWQNQYEIDIYQLGWRIGGKGASGRNAAMGQRIEEHGLHIWFGFYQNAFALMRKAYEELDRPAQAPLATFLDAFKPHHFIVLQEEINGRSVSWPI";
            const string modStr = "Acetyl 0,Oxidation 55,Dehydro 141,Dehydro 181";

            
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
            var proteinMass = sequence.Mass;

            Console.WriteLine("Mass = {0}", proteinMass);
            

            const int maxCharge = 15;
            const int minCharge = 1;
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
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

            var comparer = new FilteredProteinMassBinning(aaSet, 50000, 28);
            Console.WriteLine("{0}\t{1}",comparer.NumberOfBins, comparer.GetBinNumber(proteinMass));
            var graphFactory = new ProteinScoringGraphFactory(comparer, aaSet);

            Console.WriteLine(proteinMass);
            var stopwatch = Stopwatch.StartNew();
            var graph = graphFactory.CreateScoringGraph(productSpec, proteinMass);
            Console.WriteLine(@"graph generation elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            
            stopwatch.Restart();
            var gf = new GeneratingFunction(comparer.NumberOfBins);
            gf.ComputeGeneratingFunction(graph);
            Console.WriteLine(@"computing generation function = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            for (var score = 1; score <= gf.MaximumScore; score++)
            {
                var specEvalue = gf.GetSpectralEValue(score);
                Console.WriteLine("{0} : {1}", score, specEvalue);                                
            }


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
