using System;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestLcMsMap
    {
        public void TestRunningTime()
        {
            const string testRawFilePath = @"\\protoapps\UserData\Sangtae\Yufeng\raw\yufeng_column_test2.raw";
            var run = PbfLcMsRun.GetLcMsRun(testRawFilePath, MassSpecDataType.XCaliburRun, 0, 1.4826);
            var map = new ChargeLcScanMatrix(run);
            var scanNums = map.GetMatchingMs2ScanNums(12377.46179).ToList();
            Console.WriteLine("NumScanNums: " + scanNums.Count);
            Console.WriteLine(string.Join("\t", scanNums));
        }
    }
}
