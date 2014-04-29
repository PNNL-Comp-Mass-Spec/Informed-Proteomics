using System;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using InformedProteomics.TopDown.Execution;
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
            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_ShewIntact_2ug_3k_CID_4Apr14_Bane_PL011402.raw";
            const string dbFilePath = @"\\protoapps\UserData\Sangtae\TestData\Databases\ID_002216_235ACCEA.fasta";

            TestTopDownSearch(specFilePath, dbFilePath, true);
        }

        [Test]
        public void TestTopDownSearch(string specFilePath, string dbFilePath, bool considerInternalCleavages)
        {
            // Search parameters
            const int maxNumNTermCleavages = 30; // 30
            const int maxNumCTermCleavages = 0;
            const int minSequenceLength = 30; // 7
            const int maxSequenceLength = 250; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 30; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10; // 10
            //const int numMaxModsPerProtein = 0; // 6

            const int precursorIonTolerancePpm = 10;
            const int productIonTolerancePpm = 10;

            // Configure amino acid set
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            //var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            //var searchModifications = new List<SearchModification>
            //{
            //    //pyroGluQ,
            //    //dehydroC,
            //    //cysteinylC,
            //    //deamdN,
            //    //deamdQ,
            //    //glutathioneC,
            //    //oxM
            //};
            //var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var aaSet = new AminoAcidSet();

            IcTopDownLauncher topDownLauncher;

            if (considerInternalCleavages)
            {
                topDownLauncher = new IcTopDownLauncher(
                    specFilePath,
                    dbFilePath,
                    aaSet,
                    minSequenceLength,
                    maxSequenceLength,
                    minPrecursorIonCharge,
                    maxPrecursorIonCharge,
                    minProductIonCharge,
                    maxProductIonCharge,
                    precursorIonTolerancePpm,
                    productIonTolerancePpm,
                    true
                    );
            }
            else
            {
                // N-term cleavages up to 30 residues
                topDownLauncher = new IcTopDownLauncher(
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
                    productIonTolerancePpm
                    );
            }
            topDownLauncher.RunSearch();
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
