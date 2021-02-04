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
    internal class TestReadingProMex
    {
        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        [TestCase(8472.3354, "4325,4372")]
        [TestCase(10189.5639, "4763,4772,4773,4775,4783,4787,4789,4797,4867,4877,4885,4892,4938,4942,4943,4945")]
        public void TestReadingProMexFile(double massToFind, string expectedScanNumbers)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            var promexFilePath = Path.Combine(Utils.DEFAULT_SPEC_FILES_FOLDER, "QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.ms1ft");
            var promexFile = Utils.GetTestFile(methodName, promexFilePath);

            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName);

            Console.Write("Reading ProMex results...");
            var ms1Filter = new Ms1FtFilter(run, new Tolerance(10), promexFile.FullName);

            Console.WriteLine();

            var matchingScanNums = new SortedSet<int>();

            foreach (var item in ms1Filter.GetMatchingMs2ScanNums(massToFind))
            {
                matchingScanNums.Add(item);
            }

            var scanNumList = string.Join(",", matchingScanNums);

            Console.WriteLine("Scans with mass {0}:", massToFind);
            Console.WriteLine(scanNumList);

            var expectedScanNumList = expectedScanNumbers.Split(',');

            var matchCount = 0;
            foreach (var scanNumText in expectedScanNumList)
            {
                var scanNum = int.Parse(scanNumText);

                if (!matchingScanNums.Contains(scanNum))
                {
                    Assert.Fail("Did not find scan {0} for mass {1}", scanNum, massToFind);
                }

                matchCount++;
            }

            Assert.AreEqual(matchCount, matchingScanNums.Count, "Found extra matching scan nums vs. what was expected");
        }
    }
}
