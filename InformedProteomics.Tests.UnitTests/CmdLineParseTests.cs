using System;
using System.Reflection;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Tests.Base;
using MSPathFinderT;
using NUnit.Framework;
using PbfGen;
using ProMex;
using PRISM;

namespace InformedProteomics.Tests.UnitTests
{
    [TestFixture]
    public class CmdLineParseTests
    {
        // Ignore Spelling: Cmd, ic, tda, Cys, Acet, Dehydro, Acetyl, Deamidated

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML")]
        public void TestPbfGenParse(string filePath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;

            var args1 = new string[]
            {
                "-s", mzMLFilePath
            };

            var parser = new CommandLineParser<PbfGenInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);
            var options = result.ParsedResults;
            Assert.AreEqual(mzMLFilePath, options.SourcePath);
            var valid = options.Validate();
            Assert.IsTrue(valid, "Call to options.Validate failed");

            Console.WriteLine("Source:     " + options.SourcePath);
            Console.WriteLine("Output Dir: " + options.OutputDir);
            Console.WriteLine("Start Scan: " + options.StartScan);
            Console.WriteLine("End Scan:   " + options.EndScan);
        }

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML")]
        public void TestProMexParse(string filePath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;

            var args1 = new string[]
            {
                "-i", mzMLFilePath,
                "-minCharge", "10",
                "-maxCharge", "55",
                "-featureMap", "n",
                "-score", "y"
            };

            var parser = new CommandLineParser<ProMexInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);
            var options = result.ParsedResults;
            Assert.AreEqual(mzMLFilePath, options.InputPath);
            Assert.AreEqual(10, options.MinSearchCharge);
            Assert.AreEqual(55, options.MaxSearchCharge);
            Assert.AreEqual(false, options.FeatureMapImage);
            Assert.AreEqual(true, options.ScoreReport);
            var valid = options.Validate();
            Assert.IsTrue(valid, "Call to options.Validate failed");
            options.Display();
        }

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", @"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft", @"TEST_FOLDER\MSPathFinderT\ID_003962_71E1A1D4.fasta", @"TEST_FOLDER\Databases\Mods.txt")]
        public void TestMSPathFinderTParse(string filePath, string featurePath, string dbPath, string modsPath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;
            var fastaFile = Utils.GetTestFile(methodName, dbPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var fastaFilePath = fastaFile.FullName;
            var featureFile = Utils.GetTestFile(methodName, featurePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var featureFilePath = featureFile.FullName;
            var modsFile = Utils.GetTestFile(methodName, modsPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var modsFilePath = modsFile.FullName;

            var args1 = new string[]
            {
                "-s", mzMLFilePath,
                "-d", fastaFilePath,
                "-mod", modsFilePath,
                "-feature", featureFilePath,
                "-ic", "2",
                "-tagSearch", "0",
                "-t", "5",
                "-f", "7",
                "-tda", "1",
                "-minCharge", "10",
                "-maxCharge", "45",
                "-act", "PQD"
            };

            var parser = new CommandLineParser<TopDownInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);
            var options = result.ParsedResults;
            Assert.AreEqual(mzMLFilePath, options.SpecFilePath);
            Assert.AreEqual(fastaFilePath, options.DatabaseFilePath);
            Assert.AreEqual(modsFilePath, options.ModsFilePath);
            Assert.AreEqual(featureFilePath, options.FeatureFilePath);
            Assert.AreEqual(InternalCleavageType.MultipleInternalCleavages, options.InternalCleavageMode);
            Assert.AreEqual(false, options.TagBasedSearch);
            Assert.AreEqual(5, options.PrecursorIonTolerancePpm);
            Assert.AreEqual(7, options.ProductIonTolerancePpm);
            Assert.AreEqual(DatabaseSearchMode.Both, options.TargetDecoySearchMode);
            Assert.AreEqual(10, options.MinPrecursorIonCharge);
            Assert.AreEqual(45, options.MaxPrecursorIonCharge);
            Assert.AreEqual(ActivationMethod.PQD, options.ActivationMethod);
            var valid = options.Validate(string.Empty);
            Assert.IsTrue(valid, "Call to options.Validate failed");
            options.Display(string.Empty);
        }

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", @"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft", @"TEST_FOLDER\MSPathFinderT\ID_003962_71E1A1D4.fasta", @"TEST_FOLDER\Databases\Mods.txt")]
        public void TestMSPathFinderTParseOldSearchMode(string filePath, string featurePath, string dbPath, string modsPath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;
            var fastaFile = Utils.GetTestFile(methodName, dbPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var fastaFilePath = fastaFile.FullName;
            var featureFile = Utils.GetTestFile(methodName, featurePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var featureFilePath = featureFile.FullName;
            var modsFile = Utils.GetTestFile(methodName, modsPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var modsFilePath = modsFile.FullName;

            var args1 = new string[]
            {
                "-s", mzMLFilePath,
                "-d", fastaFilePath,
                "-mod", modsFilePath,
                "-feature", featureFilePath,
                "-m", "0",
                "-tagSearch", "0",
                "-t", "5",
                "-f", "7",
                "-tda", "1",
                "-minCharge", "10",
                "-maxCharge", "45",
                "-act", "PQD"
            };

            var parser = new CommandLineParser<TopDownInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);
            var options = result.ParsedResults;
            Assert.AreEqual(mzMLFilePath, options.SpecFilePath);
            Assert.AreEqual(fastaFilePath, options.DatabaseFilePath);
            Assert.AreEqual(modsFilePath, options.ModsFilePath);
            Assert.AreEqual(featureFilePath, options.FeatureFilePath);
            Assert.AreEqual(InternalCleavageType.MultipleInternalCleavages, options.InternalCleavageMode);
            Assert.AreEqual(false, options.TagBasedSearch);
            Assert.AreEqual(5, options.PrecursorIonTolerancePpm);
            Assert.AreEqual(7, options.ProductIonTolerancePpm);
            Assert.AreEqual(DatabaseSearchMode.Both, options.TargetDecoySearchMode);
            Assert.AreEqual(10, options.MinPrecursorIonCharge);
            Assert.AreEqual(45, options.MaxPrecursorIonCharge);
            Assert.AreEqual(ActivationMethod.PQD, options.ActivationMethod);
            var valid = options.Validate(string.Empty);
            Assert.IsTrue(valid, "Call to options.Validate failed");
            options.Display(string.Empty);
        }

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", @"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft", @"TEST_FOLDER\MSPathFinderT\ID_003962_71E1A1D4.fasta", @"TEST_FOLDER\MSPathFinderT\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage.txt")]
        public void TestMSPathFinderTParseParamFile(string filePath, string featurePath, string dbPath, string paramsPath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;
            var fastaFile = Utils.GetTestFile(methodName, dbPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var fastaFilePath = fastaFile.FullName;
            var featureFile = Utils.GetTestFile(methodName, featurePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var featureFilePath = featureFile.FullName;
            var paramFile = Utils.GetTestFile(methodName, paramsPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var paramFilePath = paramFile.FullName;

            var args1 = new string[]
            {
                "-s", mzMLFilePath,
                "-d", fastaFilePath,
                "-ParamFile", paramFilePath,
                "-feature", featureFilePath
            };

            var parser = new CommandLineParser<TopDownInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);

            var options = result.ParsedResults;

            var valid = options.Validate(paramFile.FullName);

            Assert.IsTrue(valid, "Call to options.Validate failed");
            options.Display(paramFilePath);

            Assert.AreEqual(mzMLFilePath, options.SpecFilePath);
            Assert.AreEqual(fastaFilePath, options.DatabaseFilePath);
            Assert.AreEqual(featureFilePath, options.FeatureFilePath);

            Assert.AreEqual(InternalCleavageType.SingleInternalCleavage, options.InternalCleavageMode);
            Assert.AreEqual(true, options.TagBasedSearch, "Unexpected value for TagBasedSearch");

            Assert.AreEqual(10, options.PrecursorIonTolerancePpm);
            Assert.AreEqual(10, options.ProductIonTolerancePpm);
            Assert.AreEqual(DatabaseSearchMode.Both, options.TargetDecoySearchMode);

            Assert.AreEqual(21, options.MinSequenceLength);
            Assert.AreEqual(300, options.MaxSequenceLength);

            Assert.AreEqual(2, options.MinPrecursorIonCharge);
            Assert.AreEqual(30, options.MaxPrecursorIonCharge);

            Assert.AreEqual(1, options.MinProductIonCharge);
            Assert.AreEqual(20, options.MaxProductIonCharge);

            Assert.AreEqual(3000, options.MinSequenceMass);
            Assert.AreEqual(50000, options.MaxSequenceMass);

            Assert.AreEqual(1, options.MatchesPerSpectrumToReport);
            Assert.AreEqual(ActivationMethod.Unknown, options.ActivationMethod);

            Assert.AreEqual(4, options.MaxDynamicModificationsPerSequence);

            Assert.AreEqual(3, options.Modifications.Count);

            Console.WriteLine("Validate that mod index {0} is {1}", 0, "Oxidation");
            Assert.AreEqual("Oxidation", options.Modifications[0].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[0].Location);
            Assert.AreEqual('M', options.Modifications[0].TargetResidue);
            Assert.AreEqual(Modification.Oxidation, options.Modifications[0].Modification);
            Assert.AreEqual(15.99491463, options.Modifications[0].Mass, 0.0001);

            Console.WriteLine("Validate that mod index {0} is {1}", 1, "Dehydro");
            Assert.AreEqual("Dehydro", options.Modifications[1].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[1].Location);
            Assert.AreEqual('C', options.Modifications[1].TargetResidue);
            Assert.AreEqual(Modification.Dehydro, options.Modifications[1].Modification);
            Assert.AreEqual(-1.007825035, options.Modifications[1].Mass, 0.0001);

            Console.WriteLine("Validate that mod index {0} is {1}", 2, "Acetyl");
            Assert.AreEqual("Acetyl", options.Modifications[2].Name);
            Assert.AreEqual(SequenceLocation.ProteinNTerm, options.Modifications[2].Location);
            Assert.AreEqual('*', options.Modifications[2].TargetResidue);
            Assert.AreEqual(Modification.Acetylation, options.Modifications[2].Modification);
            Assert.AreEqual(42.0105647, options.Modifications[2].Mass, 0.0001);
        }

        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", @"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft", @"TEST_FOLDER\MSPathFinderT\ID_003962_71E1A1D4.fasta", @"TEST_FOLDER\MSPathFinderT\MSPF_UVPD_MetOx_STYPhos_LysMethDiMethTriMeth_NTermAcet_Formyl_Hydroxyl_Biotin_Crotonyl_NoInternalCleavage_10ppm.txt")]
        public void TestMSPathFinderTParseParamFileSTY(string filePath, string featurePath, string dbPath, string paramsPath)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var mzMLFilePath = mzMLFile.FullName;
            var fastaFile = Utils.GetTestFile(methodName, dbPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var fastaFilePath = fastaFile.FullName;
            var featureFile = Utils.GetTestFile(methodName, featurePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));
            var featureFilePath = featureFile.FullName;
            var paramFile = Utils.GetTestFile(methodName, paramsPath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));
            var paramFilePath = paramFile.FullName;

            var args1 = new string[]
            {
                "-s", mzMLFilePath,
                "-d", fastaFilePath,
                "-ParamFile", paramFilePath,
                "-feature", featureFilePath
            };

            var parser = new CommandLineParser<TopDownInputParameters>("CmdParseTests", "1.0");
            var result = parser.ParseArgs(args1);
            Assert.True(result.Success);

            var options = result.ParsedResults;

            var valid = options.Validate(paramFile.FullName);

            Assert.IsTrue(valid, "Call to options.Validate failed");
            options.Display(paramFilePath);

            Assert.AreEqual(mzMLFilePath, options.SpecFilePath);
            Assert.AreEqual(fastaFilePath, options.DatabaseFilePath);
            Assert.AreEqual(featureFilePath, options.FeatureFilePath);

            Assert.AreEqual(InternalCleavageType.NoInternalCleavage, options.InternalCleavageMode);
            Assert.AreEqual(true, options.TagBasedSearch, "Unexpected value for TagBasedSearch");

            Assert.AreEqual(15, options.PrecursorIonTolerancePpm);
            Assert.AreEqual(14, options.ProductIonTolerancePpm);
            Assert.AreEqual(DatabaseSearchMode.Both, options.TargetDecoySearchMode);

            Assert.AreEqual(20, options.MinSequenceLength);
            Assert.AreEqual(250, options.MaxSequenceLength);

            Assert.AreEqual(3, options.MinPrecursorIonCharge);
            Assert.AreEqual(32, options.MaxPrecursorIonCharge);

            Assert.AreEqual(1, options.MinProductIonCharge);
            Assert.AreEqual(15, options.MaxProductIonCharge);

            Assert.AreEqual(3200, options.MinSequenceMass);
            Assert.AreEqual(45000, options.MaxSequenceMass);

            Assert.AreEqual(2, options.MatchesPerSpectrumToReport);
            Assert.AreEqual(ActivationMethod.UVPD, options.ActivationMethod);

            Assert.AreEqual(5, options.MaxDynamicModificationsPerSequence);

            Console.WriteLine("Mod Count: " + options.Modifications.Count);
            Assert.AreEqual(17, options.Modifications.Count);

            var modIndex = 0;
            Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Oxidation");
            Assert.AreEqual("Oxidation", options.Modifications[modIndex].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[modIndex].Location);
            Assert.AreEqual('M', options.Modifications[modIndex].TargetResidue);
            Assert.AreEqual(Modification.Oxidation, options.Modifications[modIndex].Modification);
            Assert.AreEqual(15.99491463, options.Modifications[modIndex].Mass, 0.0001);

            for (var i = 1; i <= 3; i++)
            {
                char residue;
                switch (i)
                {
                    case 1:
                        residue = 'S';
                        break;
                    case 2:
                        residue = 'T';
                        break;
                    case 3:
                        residue = 'Y';
                        break;

                    default:
                        residue = '.';
                        break;
                }

                modIndex = i;
                Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Phospho");
                Assert.AreEqual("Phospho", options.Modifications[modIndex].Name);
                Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[modIndex].Location);
                Assert.AreEqual(residue, options.Modifications[modIndex].TargetResidue);
                Assert.AreEqual(Modification.Phosphorylation, options.Modifications[modIndex].Modification);
                Assert.AreEqual(79.9663309, options.Modifications[modIndex].Mass, 0.0001);
            }

            modIndex++;
            Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Acetyl");
            Assert.AreEqual("Acetyl", options.Modifications[modIndex].Name);
            Assert.AreEqual(SequenceLocation.ProteinNTerm, options.Modifications[modIndex].Location);
            Assert.AreEqual('*', options.Modifications[modIndex].TargetResidue);
            Assert.AreEqual(Modification.Acetylation, options.Modifications[modIndex].Modification);
            Assert.AreEqual(42.0105647, options.Modifications[modIndex].Mass, 0.0001);

            modIndex++;
            Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Acetyl");
            Assert.AreEqual("Acetyl", options.Modifications[modIndex].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[modIndex].Location);
            Assert.AreEqual('K', options.Modifications[modIndex].TargetResidue);
            Assert.AreEqual(Modification.Acetylation, options.Modifications[modIndex].Modification);
            Assert.AreEqual(42.0105647, options.Modifications[modIndex].Mass, 0.0001);

            modIndex++;
            Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Methyl");
            Assert.AreEqual("Methyl", options.Modifications[modIndex].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[modIndex].Location);
            Assert.AreEqual('K', options.Modifications[modIndex].TargetResidue);
            Assert.AreEqual(Modification.Methylation, options.Modifications[modIndex].Modification);
            Assert.AreEqual(14.01565, options.Modifications[modIndex].Mass, 0.0001);

            modIndex = 16;
            Console.WriteLine("Validate that mod index {0} is {1}", modIndex, "Deamidated");
            Assert.AreEqual("Deamidated", options.Modifications[modIndex].Name);
            Assert.AreEqual(SequenceLocation.Everywhere, options.Modifications[modIndex].Location);
            Assert.AreEqual('Q', options.Modifications[modIndex].TargetResidue);
            Assert.AreEqual(Modification.Deamidation, options.Modifications[modIndex].Modification);
            Assert.AreEqual(0.984015595, options.Modifications[modIndex].Mass, 0.0001);
        }
    }
}
