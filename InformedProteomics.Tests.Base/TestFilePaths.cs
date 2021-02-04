using System;
using System.IO;
using System.Reflection;
using NUnit.Framework;

namespace InformedProteomics.Tests.Base
{
    [TestFixture]
    public class TestFilePaths
    {
        [Test]
        [TestCase(@"TEST_FOLDER\TopDown\ProductionQCShew\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw", 694084380)]
        [TestCase(@"TEST_FOLDER\TopDown\ProductionQCShew\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.pbf", 804902631)]
        [TestCase(@"TEST_FOLDER\TopDown\E_coli_iscU_60_mock.raw", 947525331)]
        [TestCase(@"TEST_FOLDER\TopDown\SBEP_STM_001_02272012_Aragon.raw", 1524969540)]
        [TestCase(@"TEST_FOLDER\TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw", 680290096)]
        [TestCase(@"TEST_FOLDER\TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.pbf", 344665109)]
        [TestCase(@"TEST_FOLDER\MSPathFinderT\ID_002216_235ACCEA.fasta", 1671864)]
        [TestCase(@"TEST_FOLDER\MSPathFinderT\ID_005133_8491EFA2.fasta", 3077075)]
        [TestCase(@"TEST_FOLDER\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", 62492003)]
        [TestCase(@"TEST_FOLDER\PNNLOmicsElementData.xml", 63982)]
        [TestCase(@"TEST_FOLDER\SpecFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt_IcTda.tsv", 76189)]
        [TestCase(@"TEST_FOLDER\TopDown\ProductionQCShew\Mods.txt", 1892)]
        [TestCase(@"TEST_FOLDER\IdFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTarget.tsv", 684649)]
        [TestCase(@"TEST_FOLDER\IdFiles\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcDecoy.tsv", 356790)]
        [TestCase(@"\\proto-2\UnitTest_Files\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw", 0)]
        [TestCase(@"\\proto-2\UnitTest_Files\E_coli_iscU_60_mock.raw", 0)]
        public void TestFiles(string fileNameOrPath, long expectedSizeBytes)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            string filePath;

            if (fileNameOrPath.StartsWith(@"TEST_FOLDER\"))
            {
                filePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, fileNameOrPath.Replace(@"TEST_FOLDER\", ""));
            }
            else
            {
                filePath = fileNameOrPath;
            }

            var testFile = Utils.GetTestFile(methodName, filePath);

            if (testFile == null)
            {
                Assert.Ignore("File not found: " + filePath);
            }

            Console.WriteLine("File found at {0}, size {1} bytes", testFile, testFile.Length);

            Assert.AreEqual(expectedSizeBytes, testFile.Length, "File size mismatch for {0}", testFile.FullName);
        }
    }
}
