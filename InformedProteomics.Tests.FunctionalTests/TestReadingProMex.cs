using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    class TestReadingProMex
    {

        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);


            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName);

            const string promexFileName = @"\\Proto-5\VOrbiETD02\2014_3\QC_Shew_Intact_26Sep14_Bane_C2Column3\MSP201508271107_Auto1226713\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            Console.Write("Reading ProMex results...");
            var ms1Filter = new Ms1FtFilter(run, new Tolerance(10), promexFileName);
            Console.WriteLine(string.Join(",", ms1Filter.GetMatchingMs2ScanNums(3016.6583)));
        }
    }
}
