using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
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

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            IMs1FeatureExtract extractor = new ChargeLcScanMatrix(run, minScanCharge, maxScanCharge);
            var outputFilePath = extractor.GetFeatureFile(specFilePath, minScanMass, maxScanMass);
        }
    }
}
