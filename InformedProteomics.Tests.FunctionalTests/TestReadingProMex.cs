using System;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests.FunctionalTests
{
    [TestFixture]
    class TestReadingProMex
    {
        [Test]
        public void TestReadingProMexFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_214\2014_3\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const string promexFileName = @"\\Proto-5\VOrbiETD02\2014_3\QC_Shew_Intact_26Sep14_Bane_C2Column3\MSP201508271107_Auto1226713\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            Console.Write("Reading ProMex results...");
            var ms1Filter = new Ms1FtFilter(run, new Tolerance(10), promexFileName);
            Console.WriteLine(string.Join(",", ms1Filter.GetMatchingMs2ScanNums(3016.6583)));
        }
    }
}
