using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;


namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestProMex
    {
        [Test]
        public void TestGeneratingMs1FeatureFile()
        {
            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";
            
            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 3100;
            const int maxThreads = 10;

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            IMs1FeatureExtract extractor = new ChargeLcScanMatrix(run, minScanCharge, maxScanCharge, maxThreads);
            var outputFilePath = extractor.GetFeatureFile(specFilePath, minScanMass, maxScanMass);
        }

        [Test]
        public void TestProMexFilter()
        {
            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            const string ms1FtPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath, 0.15);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }
    }
}
