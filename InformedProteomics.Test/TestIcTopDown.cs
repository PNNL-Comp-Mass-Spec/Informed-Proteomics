using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcTopDown
    {
        [Test]
        [TestCase(@"TEST_FOLDER\MSPathFinderT\ID_004530_B63BD900.fasta", 9146396, 10562511, 11333064)]
        public void TestForManyMods(string dbFile, int expectedMaxLen300, int expectedMaxLen400, int expectedMaxLen500)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var fastaFile = Utils.GetTestFile(methodName, dbFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            var indexedDb = new IndexedDatabase(new FastaDatabase(fastaFile.FullName));
            indexedDb.Read();

            var totalPeptidesMaxLen300 = indexedDb.EstimateTotalPeptides(1, 21, 300, 1, 0);
            Console.WriteLine("{0,12:N0} estimated peptides for max length 300", totalPeptidesMaxLen300);

            var totalPeptidesMaxLen400 = indexedDb.EstimateTotalPeptides(1, 21, 400, 1, 0);
            Console.WriteLine("{0,12:N0} estimated peptides for max length 400", totalPeptidesMaxLen400);

            var totalPeptidesMaxLen500 = indexedDb.EstimateTotalPeptides(1, 21, 500, 1, 0);
            Console.WriteLine("{0,12:N0} estimated peptides for max length 500", totalPeptidesMaxLen500);

            Assert.AreEqual(expectedMaxLen300, totalPeptidesMaxLen300);
            Assert.AreEqual(expectedMaxLen400, totalPeptidesMaxLen400);
            Assert.AreEqual(expectedMaxLen500, totalPeptidesMaxLen500);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestForVlad()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"D:\Research\Data\Vlad\raw\Alz_RA_C1_HCD_11012013_SW_03Nov2013.raw";
            const string dbFilePath = @"D:\Research\Data\Vlad\database\ID_004221_1C042A1F.fasta";
            //const string dbFilePath = @"D:\Research\Data\Vlad\database\HBA_MOUSE.fasta";
            const string outputDir = @"D:\Research\Data\Vlad\Ic\POPSICLETest_M1";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // ReSharper disable IdentifierTypo
            // ReSharper disable CommentTypo
            // ReSharper disable UnusedVariable

            // Configure amino acid set
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var thrToAla = new SearchModification(Modification.ThrToAla, 'T', SequenceLocation.Everywhere, false);
            var dethiomethylM = new SearchModification(Modification.Dethiomethyl, 'M', SequenceLocation.Everywhere, false);
            var deamidatedN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            var deamidatedQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);
            var serToAsn = new SearchModification(Modification.SerToAsn, 'S', SequenceLocation.Everywhere, false);
            var pyroCarbamidomethylC = new SearchModification(Modification.PyroCarbamidomethyl, 'C',
                SequenceLocation.ProteinNTerm, false);
            var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);

            // ReSharper restore UnusedVariable

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                //                glutathioneC,
                oxM,
//                dethiomethylM,
                acetylN,
                phosphoS,
                phosphoT,
                phosphoY
//                thrToAla,
//                serToAsn,
//                deamidatedN,
//                deamidatedQ,
//                pyroCarbamidomethylC
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            // ReSharper restore CommentTypo
            // ReSharper restore IdentifierTypo

            const InternalCleavageType searchMode = InternalCleavageType.SingleInternalCleavage;   // 0: all sub-sequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Target;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, searchMode);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestForYufeng()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew
            const string specFilePath = @"H:\Research\Yufeng\TopDownYufeng\raw\yufeng_column_test2.raw";
            //const string dbFilePath = @"H:\Research\Yufeng\TopDownYufeng\database\ID_002216_235ACCEA.fasta";
            const string dbFilePath = @"H:\Research\Yufeng\TopDownYufeng\database\SO_3942_Truncated.fasta";
            const string outputDir = @"H:\Research\Yufeng\TopDownYufeng\Debug";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // ReSharper disable IdentifierTypo
            // ReSharper disable CommentTypo

            // Configure amino acid set
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            ////            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            //var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            //var pyroGluQ = new SearchModification(Modification.PTeyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            //const int numMaxModsPerProtein = 0;
            //var searchModifications = new List<SearchModification>
            //{
            //    dehydroC,
            //    oxM,
            //    acetylN
            //};
            //var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            // ReSharper restore CommentTypo
            // ReSharper restore IdentifierTypo

            var aaSet = new AminoAcidSet();
            const InternalCleavageType searchMode = InternalCleavageType.NoInternalCleavage;   // 0: all sub-sequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Target;   // true: target & decoy, false: target, null: decoy

            const int minSequenceLength = 21; // 7
            const int maxSequenceLength = 500; // 1000
            const int minPrecursorIonCharge = 2; // 3
            const int maxPrecursorIonCharge = 50; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 20; // 15
            const double minSequenceMass = 3000.0;
            const double maxSequenceMass = 50000.0;

            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet,
                minSequenceLength, maxSequenceLength,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge,
                minSequenceMass, maxSequenceMass,
                tda, searchMode
                );
        }

        [Test]
        [Category("PNL_Domain")]
        [Category("Local_Testing")]
        public void TestForSystemsBiologyData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //// Salmonella
            var specFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\SBEP_STM_001_02272012_Aragon.raw");
            var dbFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"MSPathFinderT\ID_002166_F86E3B2F.fasta");
            var outputDir = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"Results\Mod_M2");

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            if (!Directory.Exists(outputDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, outputDir);
            }

            // ReSharper disable IdentifierTypo

            // Configure amino acid set
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);

            // ReSharper restore IdentifierTypo

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                oxM,
                acetylN,
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const InternalCleavageType searchMode = InternalCleavageType.NoInternalCleavage;   // 0: all sub-sequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Both;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, searchMode);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestForAaronData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownAaron\raw\MTB_intact_1.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownAaron\database\ID_003121_998584F8.fasta";
            const string outputDir = @"C:\cygwin\home\kims336\Data\TopDownAaron\Ic\Mode1_07";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // Configure amino acid set
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var tevFp2C = new SearchModification(Modification.TevFp2, 'S', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                //glutathioneC,
                //nitrosylC,
                //nethylmaleimideC,
                oxM,
                acetylN,
                tevFp2C
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            const InternalCleavageType searchMode = InternalCleavageType.SingleInternalCleavage;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Both;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, searchMode);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestForJiaData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew
            //const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            // var dbFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"MSPathFinderT\ID_002216_235ACCEA.fasta");
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\Test.fasta";

            // Jia's data
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\ID_003962_71E1A1D4.fasta";
            const string outputDir = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\D1_1_Mode1";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // Configure amino acid set
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);
            //var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                nitrosylC,
                nethylmaleimideC,
                oxM,
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const InternalCleavageType searchMode = InternalCleavageType.SingleInternalCleavage;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Both;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, searchMode);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestForQcShew()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew
            const string specFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string dbFilePath = @"D:\MSPathFinder\Fasta\ID_002216_235ACCEA.fasta";
            const string outputDir = @"D:\MassSpecFiles\training\test";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
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

            const InternalCleavageType searchMode = InternalCleavageType.NoInternalCleavage;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            const DatabaseSearchMode tda = DatabaseSearchMode.Both;   // true: target & decoy, false: target, null: decoy
            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, searchMode);
        }

        private void TestTopDownSearch(string specFilePath, string dbFilePath, string outputDir, AminoAcidSet aaSet,
            DatabaseSearchMode tda, InternalCleavageType searchMode)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int minSequenceLength = 21; // 7
            const int maxSequenceLength = 500; // 1000
            const int minPrecursorIonCharge = 2; // 3
            const int maxPrecursorIonCharge = 60; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 20; // 15
            const double minSequenceMass = 3000.0;
            const double maxSequenceMass = 50000.0;

            TestTopDownSearch(specFilePath, dbFilePath, outputDir, aaSet,
                minSequenceLength, maxSequenceLength,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge,
                minSequenceMass, maxSequenceMass,
                tda, searchMode
                );
        }

        private void TestTopDownSearch(string specFilePath, string dbFilePath, string outputDir, AminoAcidSet aaSet,
            int minSequenceLength, int maxSequenceLength,
            int minPrecursorIonCharge, int maxPrecursorIonCharge,
            int minProductIonCharge, int maxProductIonCharge,
            double minSequenceMass, double maxSequenceMass,
            DatabaseSearchMode tda, InternalCleavageType searchMode)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // Search parameters
            const int maxNumNTermCleavages = 1; // 30
            const int maxNumCTermCleavages = 0;
            const int precursorIonTolerancePpm = 10;
            const int productIonTolerancePpm = 10;

            var topDownOptions = new MsPfParameters(
                    specFilePath,
                    dbFilePath,
                    outputDir,
                    aaSet, "")
                    {
                    MinSequenceLength = minSequenceLength,
                    MaxSequenceLength = maxSequenceLength,
                    MaxNumNTermCleavages = maxNumNTermCleavages,
                    MaxNumCTermCleavages = maxNumCTermCleavages,
                    MinPrecursorIonCharge = minPrecursorIonCharge,
                    MaxPrecursorIonCharge = maxPrecursorIonCharge,
                    MinProductIonCharge = minProductIonCharge,
                    MaxProductIonCharge = maxProductIonCharge,
                    MinSequenceMass = minSequenceMass,
                    MaxSequenceMass = maxSequenceMass,
                    PrecursorIonTolerancePpm = precursorIonTolerancePpm,
                    ProductIonTolerancePpm = productIonTolerancePpm,
                    TargetDecoySearchMode = tda,
                    InternalCleavageMode = searchMode,
                    };

            var topDownLauncher = new IcTopDownLauncher(topDownOptions);
            //topDownLauncher.ForceParallel = true;
            //topDownLauncher.MaxNumThreads = -1;

            topDownLauncher.RunSearch(0.7);
            //topDownLauncher.RunIntactProteinSearch();
        }

        [Test]
        [Category("Local_Testing")]
        public void TestPrSm()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownYufeng\raw\yufeng_column_test2.raw";
            //const string annotation =
            //    "_.MKTKLSVLSAAMLAATLTMMPAVSQAAIPQSVEGQSIPSLAPMLERTTPAVVSVAVSGTHVSKQRVPDVFRYFFGPNAPQEQVQERPFRGLGSGVIIDADKGYIVTNNHVIDGADDIQVG" +
            //    "LHDGREVKAKLIGTDSESDIALLQIEAKNLVAIKTSDSDELRVGDFAVAIGNPFGLGQTV" +
            //    "TSGIVSALGRSGLGIEMLENFIQTDAAINSGNSGGALVNLKGELIGINTAIVAPNGGNVG" +
            //    "IGFAIPANMVKNLIAQIAEHGEVRRGVLGIAGRDLDSQLAQGFGLDTQHGGFVNEVSAGS" +
            //    "AAEKAGIKAGDIIVSVDGRAIKSFQELRAKVATMGAGAKVELGLIRDGDKKTVNVTLGEA" +
            //    "NQTTEKAAGAVHPMLQGASLENASKGVEITDVAQGSPAAMSGLQKGDLIVGINRTAVKDL" +
            //    "KSLKELLKDQEGAVALKIVRGKSMLYLVLR._";
            //var aaSet = new AminoAcidSet();

            //const int charge = 60;
            //const int ms2ScanNum = 46661;

            const string specFilePath = @"D:\Research\Data\Jon\AH_SF_mouseliver_3-1_Intact_2_6Feb14_Bane_PL011402.raw";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const int ms2ScanNum = 19011;
            const int charge = 7;
            const string annotation = "_.SKVSFKITLTSDPRLPYKVLSVPESTPFTAVLKFAAEEFKVPAATSAIITNDGIGINPAQTAGNVFLKHGSELRIIPRDRVGSC._";

            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, true);
            var modVal = Modification.RegisterAndGetModification("AddVal", new Composition(5, 9, 1, 1, 0));
            var searchMods = AminoAcid.StandardAminoAcidCharacters.Select(residue => new SearchModification(modVal, residue, SequenceLocation.Everywhere, false)).ToList();
            searchMods.Add(acetylN);
            const int numMaxModsPerProtein = 1;
            var aaSet = new AminoAcidSet(searchMods, numMaxModsPerProtein);

            var graph = SequenceGraph.CreateGraph(aaSet, annotation);
            Console.WriteLine("NumProteoforms: " + graph.GetNumProteoformCompositions());

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, 1.4826, 1.4826);
            var ms2Scorer = new ProductScorerBasedOnDeconvolutedSpectra(run, 1, 15);
            ms2Scorer.GetScorer(ms2ScanNum);
            var scorer = ms2Scorer.GetMs2Scorer(ms2ScanNum);
            Assert.NotNull(scorer, "Scorer is null!");

            for (var i = 0; i < graph.GetNumProteoformCompositions(); i++)
            {
                graph.SetSink(i);
                Console.WriteLine("ModComb: " + graph.GetModificationCombinations()[i]);
                var score = graph.GetFragmentScore(scorer);
                Console.WriteLine("Fast search score: " + score);
                var composition = graph.GetSinkSequenceCompositionWithH2O();

                var informedScorer = new InformedTopDownScorer(run, aaSet, 1, 30, new Tolerance(10));
                var refinedScore = informedScorer.GetScores(AminoAcid.ProteinNTerm, SimpleStringProcessing.GetStringBetweenDots(annotation), AminoAcid.ProteinCTerm, composition, charge, ms2ScanNum);
                Console.WriteLine("Modifications: {0}", refinedScore.Modifications);
                Console.WriteLine("Composition: {0}", composition);
                Console.WriteLine("RefinedScores: {0}", refinedScore);
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void TestMsAlignRescoring()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var specFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw");
            const string msAlignResultPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\NoMod.tsv";
            const string outputPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\NoMod_Rescored.tsv";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(msAlignResultPath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, msAlignResultPath);
            }

            var tolerance = new Tolerance(10.0);

            var rescorer = new MsAlignRescorer(specFilePath, msAlignResultPath, outputPath, tolerance, 0.7, 1, 15);
            Console.WriteLine("Done");
        }

        [Test]
        [Category("Local_Testing")]
        public void TestIcRescoring()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // var specFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, @"QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw");
            // var icResultPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402_Map07_Re.icdresult";
            // var outputPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\raw\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402_Map07_Re_Rescored.icdresult";

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            const string icResultPath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.ictresult";
            const string outputPath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1_Rescored.ictresult";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(icResultPath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, icResultPath);
            }

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

            var reScorer = new IcRescorer(specFilePath, icResultPath, outputPath, aaSet, tolerance, 0.7);
            Console.WriteLine("Done");
        }

        //[Test]
        //public void TestHistoneSearch()
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
