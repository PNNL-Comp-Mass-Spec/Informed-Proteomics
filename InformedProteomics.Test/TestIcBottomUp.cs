using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.BottomUp.Execution;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestIcBottomUp
    {
        [Test]
        public void TestChaoChaoWhim()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFileDir = @"D:\Research\Data\ChaoChao\WHIM\raw";
            if (!Directory.Exists(specFileDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, specFileDir);
            }

            foreach (var specFilePath in Directory.GetFiles(specFileDir, "*.raw"))
            {
                Console.WriteLine("****** Processing {0}", specFilePath);
                TestChaoChao(specFilePath);
            }
        }

        public void TestChaoChao(string specFilePath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const string dbFilePath = @"D:\Research\Data\ChaoChao\database\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            const string outputDir = @"D:\Research\Data\ChaoChao\Ic\";

            // Configure amino acid set
            //var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 0;
            var searchModifications = new List<SearchModification>
            {
                //carbamidomethylC,
                //acetylN,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            bool? tda = true;   // true: target & decoy, false: target, null: decoy

            const int minSequenceLength = 7; // 7
            const int maxSequenceLength = 150; // 1000
            const int minPrecursorIonCharge = 1; // 3
            const int maxPrecursorIonCharge = 30; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 15; // 15
            const double precursorIonTolerancePpm = 10;
            const double productIonTolerancePpm = 10;
            const double corrThreshold = 0.7;

            var bottomUpLauncher = new IcBottomUpLauncher(
                specFilePath,
                dbFilePath,
                outputDir,
                aaSet,
                null)
            {
                MinSequenceLength = minSequenceLength,
                MaxSequenceLength = maxSequenceLength,
                MinPrecursorIonCharge = minPrecursorIonCharge,
                MaxPrecursorIonCharge = maxPrecursorIonCharge,
                MinProductIonCharge = minProductIonCharge,
                MaxProductIonCharge = maxProductIonCharge,
                PrecursorIonTolerancePpm = precursorIonTolerancePpm,
                ProductIonTolerancePpm = productIonTolerancePpm,
                RunTargetDecoyAnalysisBool = tda,
                NumTolerableTermini = 0
            };

            bottomUpLauncher.RunSearch(corrThreshold);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestEdrn()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFileDir = @"H:\Research\EDRN\RawFiles\DIA";
            if (!Directory.Exists(specFileDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, specFileDir);
            }

            foreach (var specFilePath in Directory.GetFiles(specFileDir, "*DIA*.raw"))
            {
                Console.WriteLine("****** Processing {0}", specFilePath);
                TestEdrn(specFilePath);
            }
        }

        private void TestEdrn(string specFilePath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const string dbFilePath = @"H:\Research\EDRN\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            const string outputDir = @"H:\Research\EDRN\Ic";

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // Configure amino acid set
            var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 0;
            var searchModifications = new List<SearchModification>
            {
                carbamidomethylC,
                //acetylN,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const int ntt = 2;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            bool? tda = true;   // true: target & decoy, false: target, null: decoy
            TestBottomUpSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, ntt);
        }

        [Test]
        public void TestMaccossDia()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFileDir = @"D:\Research\Data\UW\QExactive\";
            if (!Directory.Exists(specFileDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, specFileDir);
            }

            foreach (var specFilePath in Directory.GetFiles(specFileDir, "*DIA*.raw"))
            {
                Console.WriteLine("****** Processing {0}", specFilePath);
                TestMaccoss(specFilePath);
            }
        }

        [Test]
        public void TestMaccossDda()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFileDir = @"D:\Research\Data\UW\QExactive\";
            if (!Directory.Exists(specFileDir))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, specFileDir);
            }

            foreach (var specFilePath in Directory.GetFiles(specFileDir, "*DDA*.raw"))
            {
                Console.WriteLine("****** Processing {0}", specFilePath);
                TestMaccoss(specFilePath);
            }
        }

        public void TestMaccoss(string specFilePath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const string dbFilePath = @"D:\Research\Data\UW\QExactive\M_musculus_Uniprot_withContam.fasta";
            const string outputDir = @"D:\Research\Data\UW\QExactive\Ic_NTT1_Mass";

            if (!File.Exists(dbFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dbFilePath);
            }

            // Configure amino acid set
            var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 0;
            var searchModifications = new List<SearchModification>
            {
                carbamidomethylC,
                //acetylN,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const int ntt = 1;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            bool? tda = true;   // true: target & decoy, false: target, null: decoy
            TestBottomUpSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, ntt);
        }

        [Test]
        public void TestQcShewQExactive()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew QE
            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\ID_003456_9B916A8B.fasta";
            const string outputDir = @"C:\cygwin\home\kims336\Data\QCShewQE\Ic_NTT2_Mass";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            // Configure amino acid set
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);

            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 0;
            var searchModifications = new List<SearchModification>
            {
                //carbamidomethylC,
                acetylN,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const int ntt = 2;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            bool? tda = true;   // true: target & decoy, false: target, null: decoy
            TestBottomUpSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, ntt);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestDdaPlus()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew QE
            const string specFilePath = @"H:\Research\DDAPlus\raw\20140701_yeast_DDA_01.raw";
            const string dbFilePath = @"H:\Research\DDAPlus\database\Yeast_SGD_withContam.fasta";
            const string outputDir = @"H:\Research\DDAPlus\Test";

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
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);

            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            const int numMaxModsPerPeptide = 2;
            var searchModifications = new List<SearchModification>
            {
                //carbamidomethylC,
                acetylN,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerPeptide);

            const int ntt = 2;   // 0: all subsequences, 1: close to N- or C-term, 2: close to N- and C-term
            bool? tda = true;   // true: target & decoy, false: target, null: decoy
            TestBottomUpSearch(specFilePath, dbFilePath, outputDir, aaSet, tda, ntt);
        }

        public void TestBottomUpSearch(string specFilePath, string dbFilePath, string outputDir, AminoAcidSet aaSet, bool? tda, int ntt, double corrThreshold = 0.3)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // Search parameters
            const int minSequenceLength = 6; // 7
            const int maxSequenceLength = 40; // 1000
            const int minPrecursorIonCharge = 1; // 3
            const int maxPrecursorIonCharge = 4; // 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 2; // 15
            const int precursorIonTolerancePpm = 10;
            const int productIonTolerancePpm = 10;

            var enzyme = Enzyme.Trypsin;

            var bottomUpLauncher = new IcBottomUpLauncher(
                specFilePath,
                dbFilePath,
                outputDir,
                aaSet,
                enzyme)
            {
                MinSequenceLength = minSequenceLength,
                MaxSequenceLength = maxSequenceLength,
                MinPrecursorIonCharge = minPrecursorIonCharge,
                MaxPrecursorIonCharge = maxPrecursorIonCharge,
                MinProductIonCharge = minProductIonCharge,
                MaxProductIonCharge = maxProductIonCharge,
                PrecursorIonTolerancePpm = precursorIonTolerancePpm,
                ProductIonTolerancePpm = productIonTolerancePpm,
                RunTargetDecoyAnalysisBool = tda,
                NumTolerableTermini = ntt
            };

            bottomUpLauncher.RunSearch(corrThreshold);
            //topDownLauncher.RunIntactProteinSearch();
        }
    }
}
