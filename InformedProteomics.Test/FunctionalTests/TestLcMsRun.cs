using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestLcMsRun
    {
        public const string TestRawFilePath = @"\\protoapps\UserData\Sangtae\TestData\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";
        public const string TestPbfFilePath = @"\\protoapps\UserData\Sangtae\TestData\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.pbf";
        public const string TestTopDownRawFilePathEtd = @"\\protoapps\UserData\Sangtae\TestData\E_coli_iscU_60_mock.raw";
        public const string TestTopDownRawFilePathCid = @"\\protoapps\UserData\Sangtae\TestData\SBEP_STM_001_02272012_Aragon.raw";

        //"\\protoapps\UserData\Sangtae\TopDownQCShew\raw";
        public void TestReadingScanNums()
        {
            var run = LcMsRun.GetLcMsRun(TestRawFilePath, MassSpecDataType.XCaliburRun);

            var msLevel = new Dictionary<int, int>();

            for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            {
                msLevel[scanNum] = run.GetMsLevel(scanNum);
            }

            for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            {
                var spec = run.GetSpectrum(scanNum);
                Assert.True(spec.MsLevel == msLevel[scanNum]);

                if (spec.MsLevel == 2)
                {
                    var precursorScanNum = 0;
                    for (var prevScan = scanNum - 1; prevScan >= run.MinLcScan; prevScan--)
                    {
                        if (run.GetMsLevel(prevScan) == 1)
                        {
                            precursorScanNum = prevScan;
                            break;
                        }
                    }
                    Assert.True(run.GetPrecursorScanNum(scanNum) == precursorScanNum);

                    var nextScanNum = run.MaxLcScan+1;
                    for (var nextScan = scanNum + 1; nextScan <= run.MaxLcScan; nextScan++)
                    {
                        if (run.GetMsLevel(nextScan) == 1)
                        {
                            nextScanNum = nextScan;
                            break;
                        }
                    }
                    if (run.GetNextScanNum(scanNum) != nextScanNum)
                    {
                        Console.WriteLine("{0}\t{1}\t{2}", scanNum, run.GetNextScanNum(scanNum), nextScanNum);
                    }
                    Assert.True(run.GetNextScanNum(scanNum) == nextScanNum);
                }
            }
            Assert.True(run.GetNextScanNum(31151) == 31153);
            Console.WriteLine(run.GetNextScanNum(89));
        }

        //[Test]
        //public void TestReadingPbf()
        //{
        //    var sw = new System.Diagnostics.Stopwatch();

        //    sw.Start();
        //    var run = new LcMsRun(new PbfReader(TestPbfFilePath));
        //    Console.WriteLine(run.MaxLcScan);
        //    var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
        //    Console.WriteLine(@"Done. {0:f4} sec", sec);
        //}

        [Test]
        public void TestParsingSpectrumFile()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            var run = LcMsRun.GetLcMsRun(TestTopDownRawFilePathEtd, MassSpecDataType.XCaliburRun);

            const int scanNum = 810;
            var spec = run.GetSpectrum(scanNum) as ProductSpectrum;
            if (spec != null)
            {
                spec.Display();
                var precursorInfo = spec.IsolationWindow;
                Console.WriteLine("ActivationMethod: {0}", spec.ActivationMethod);
                Console.WriteLine("PrecursorScan: {0}", run.GetPrecursorScanNum(spec.ScanNum));
                Console.WriteLine("IsolationWindowTargetMz: {0}", precursorInfo.IsolationWindowTargetMz);
                Console.WriteLine("IsolationWindowLowerOffset: {0}", precursorInfo.IsolationWindowLowerOffset);
                Console.WriteLine("IsolationWindowUpperOffset: {0}", precursorInfo.IsolationWindowUpperOffset);
                Console.WriteLine("MsLevel: {0}", run.GetMsLevel(scanNum));
            }

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestNoiseFiltration()
        {
            
        }
    }
}
