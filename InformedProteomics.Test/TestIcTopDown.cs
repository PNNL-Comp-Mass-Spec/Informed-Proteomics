using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcTopDown
    {
        [Test]
        public void TestTopDownSearch()
        {
            //// Salmonella
            //const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\SBEP_STM_001_02272012_Aragon.raw";
            //const string dbFilePath = @"\\protoapps\UserData\Sangtae\TestData\Databases\ID_002166_F86E3B2F.fasta";

            // QC_Shew
            //const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            //const string dbFilePath = @"\\protoapps\UserData\Sangtae\TestData\Databases\ID_002216_235ACCEA.fasta";
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\Test.fasta";

            // var aaSet = new AminoAcidSet();

            // Jia's data
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\ID_003962_71E1A1D4.fasta";

            // Configure amino acid set
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                nitrosylC,
                nethylmaleimideC,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const int searchMode = 1;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            bool? tda = true;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, aaSet, tda, searchMode);
        }

        [Test]
        public void TestTopDownSearch(string specFilePath, string dbFilePath, AminoAcidSet aaSet, bool? tda, int searchMode)
        {
            // Search parameters
            const int maxNumNTermCleavages = 1; // 30
            const int maxNumCTermCleavages = 0;
            const int minSequenceLength = 21; // 7
            const int maxSequenceLength = 300; // 1000
            const int minPrecursorIonCharge = 2; // 3
            const int maxPrecursorIonCharge = 30; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 15; // 15

            const int precursorIonTolerancePpm = 10;
            const int productIonTolerancePpm = 10;

            var topDownLauncher = new IcTopDownLauncher(
                    specFilePath,
                    dbFilePath,
                    aaSet,
                    minSequenceLength,
                    maxSequenceLength,
                    maxNumNTermCleavages,
                    maxNumCTermCleavages,
                    minPrecursorIonCharge,
                    maxPrecursorIonCharge,
                    minProductIonCharge,
                    maxProductIonCharge,
                    precursorIonTolerancePpm,
                    productIonTolerancePpm,
                    tda,
                    searchMode
                    );

            topDownLauncher.RunSearch(0.7);
            //topDownLauncher.RunIntactProteinSearch();
        }

        [Test]
        public void TestPrSm()
        {
            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            const string annotation =
                "M.AAMTVQLDIVSAESSIFSGRVASLQVTGSEGELGIMHGHAPLLSYIKPGMARIVKQDGNEEVFYLSGGLLEVQPSSVSVLADVVMRAKDIDEQAALEAKRRAEAHMATAGADFNYDAAMVELAKAMAQLRVVETIKKNIAR._";
            const int charge = 15;
            const int ms2ScanNum = 5595;

            var aaSet = new AminoAcidSet();

            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var ms2Scorer = new ProductScorerBasedOnDeconvolutedSpectra(run, 1, 15);
            ms2Scorer.DeconvoluteProductSpectra();
            var scorer = ms2Scorer.GetMs2Scorer(ms2ScanNum);

            var graph = SequenceGraph.CreateGraph(aaSet, annotation);
            graph.SetSink(0, 0);
            var score = graph.GetScore(charge, scorer);
            Console.WriteLine("Fast search score: " + score);
            var composition = graph.GetSinkSequenceCompositionWithH2O();

            var informedScorer = new InformedScorer(run, new AminoAcidSet(), 1, 15, new Tolerance(10));
            var refinedScore = informedScorer.GetScores(SimpleStringProcessing.GetStringBetweenDots(annotation), composition, charge, ms2ScanNum).Ms2Score;
            Console.WriteLine("RefinedScore: {0}", refinedScore);
        }

        [Test]
        public void TestMsAlignRescoring()
        {
            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            const string msAlignResultPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\NoMod.tsv";
            const string outputPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\NoMod_Rescored.tsv";
            var tolerance = new Tolerance(10.0);

            var rescorer = new MsAlignRescorer(specFilePath, msAlignResultPath, outputPath, tolerance, 0.7, 1, 15);
            Console.WriteLine("Done");
        }

        [Test]
        public void TestIcRescoring()
        {
            //const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            //const string icResultPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402_Map07_Re.icdresult";
            //const string outputPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402_Map07_Re_Rescored.icdresult";

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            const string icResultPath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.icdresult";
            const string outputPath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1_Rescored.icdresult";
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                nitrosylC,
                nethylmaleimideC,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            var tolerance = new Tolerance(10.0);

            var rescorer = new IcRescorer(specFilePath, icResultPath, outputPath, aaSet, tolerance, 0.7);
            Console.WriteLine("Done");
        }

        [Test]
        public void ExtractProteinSequences()
        {
            const string fastaFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\ID_003962_71E1A1D4.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            fastaDb.Read();

            const string proteinFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\Proteins.txt";
            var proteins = File.ReadAllLines(proteinFilePath);
            foreach (var protein in proteins)
            {
                var token = protein.Split();
                if (token.Length < 1) continue;
                var proteinId = protein.Split()[0];
                var proteinSequence = fastaDb.GetProteinSequence(proteinId);
                Assert.IsTrue(proteinSequence != null);
                Console.WriteLine(">" + protein);
                Console.WriteLine(proteinSequence);
            }
        }

        //[Test]
        //public void TestHistonSearch()
        //{
        //    const bool isDecoy = true;
        //    // Search parameters
        //    const int maxNumNTermCleavages = 1;
        //    const int maxNumCTermCleavages = 0;
        //    const int minLength = 7; // 7
        //    const int maxLength = 1000; // 1000
        //    const int minPrecursorIonCharge = 3; // 3
        //    const int maxPrecursorIonCharge = 40; // 67
        //    const int minProductIonCharge = 1; // 1
        //    const int maxProductIonCharge = 10; // 10
        //    const int numMaxModsPerProtein = 11; // 6

        //    var precursorTolerance = new Tolerance(15);
        //    var productIonTolerance = new Tolerance(15);

        //    const string dbFilePath = @"D:\Research\Data\TopDownHistone\HistoneH4.fasta";

        //    //            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
        //    const string specFilePath = @"D:\Research\Data\TopDownHistone\071210_070610His0Gy070210H4_H061010A.raw";

        //    var acetylR = new SearchModification(Modification.Acetylation, 'R', SequenceLocation.Everywhere, false);
        //    var acetylK = new SearchModification(Modification.Acetylation, 'K', SequenceLocation.Everywhere, false);
        //    var methylR = new SearchModification(Modification.Methylation, 'R', SequenceLocation.Everywhere, false);
        //    var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
        //    var diMethylR = new SearchModification(Modification.DiMethylation, 'R', SequenceLocation.Everywhere, false);
        //    var diMethylK = new SearchModification(Modification.DiMethylation, 'K', SequenceLocation.Everywhere, false);
        //    var triMethylR = new SearchModification(Modification.TriMethylation, 'R', SequenceLocation.Everywhere, false);
        //    var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
        //    var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
        //    var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);

        //    var searchModifications = new List<SearchModification>
        //    {
        //        acetylR,
        //        acetylK,
        //        methylR,
        //        methylK,
        //        diMethylR,
        //        diMethylK,
        //        triMethylR,
        //        phosphoS,
        //        phosphoT,
        //        phosphoY
        //    };
        //    var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
        //}
    }
}
