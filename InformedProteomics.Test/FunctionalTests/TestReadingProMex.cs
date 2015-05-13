using System;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    class TestReadingProMex
    {
        [Test]
        public void TestReadingProMexFile()
        {
            const string rawFilePath = @"H:\Research\GlycoTopDown\raw\User_sample_test_MWCO_02262016.raw";
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const string promexFileName = @"H:\Research\GlycoTopDown\raw\User_sample_test_MWCO_02262016.ms1ft";
            Console.Write("Reading ProMex results...");
            var ms1Filter = new Ms1FtFilter(run, new Tolerance(10), promexFileName, 0.15);
            Console.WriteLine(string.Join(",", ms1Filter.GetMatchingMs2ScanNums(10000.0)));
        }
    }
}
