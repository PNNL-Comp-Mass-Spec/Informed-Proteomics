using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
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
            Assert.IsTrue(valid);
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
            Assert.IsTrue(valid);
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
            var valid = options.Validate();
            Assert.IsTrue(valid);
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
            var valid = options.Validate();
            Assert.IsTrue(valid);
        }
    }
}
