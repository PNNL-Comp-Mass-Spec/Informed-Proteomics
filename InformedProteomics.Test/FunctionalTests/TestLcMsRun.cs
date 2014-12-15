using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestLcMsRun
    {
        public const string TestRawFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";
        public const string TestPbfFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.pbf";
        public const string TestTopDownRawFilePathEtd = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\E_coli_iscU_60_mock.raw";
        public const string TestTopDownRawFilePathCid = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\SBEP_STM_001_02272012_Aragon.raw";
        public const string TestQExactiveRawFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
        public const string TestQExactivePbfFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.pbf";

        //"\\protoapps\UserData\Sangtae\TopDownQCShew\raw";
        public void TestReadingScanNums()
        {
            var run = InMemoryLcMsRun.GetLcMsRun(TestRawFilePath, MassSpecDataType.XCaliburRun);

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

        [Test]
        public void TestParsingSpectrumFile()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            var run = InMemoryLcMsRun.GetLcMsRun(TestTopDownRawFilePathCid, MassSpecDataType.XCaliburRun);

            const int scanNum = 425;
            var spec = run.GetSpectrum(scanNum) as ProductSpectrum;
            if (spec != null)
            {
                spec.Display();
                var precursorInfo = spec.IsolationWindow;
                Console.WriteLine("ActivationMethod: {0}", spec.ActivationMethod);
                Console.WriteLine("Rt: {0}", spec.ElutionTime);
                Console.WriteLine("PrecursorScan: {0}", run.GetPrecursorScanNum(spec.ScanNum));
                Console.WriteLine("IsolationWindowTargetMz: {0}", precursorInfo.IsolationWindowTargetMz);
                Console.WriteLine("IsolationWindowLowerOffset: {0}", precursorInfo.IsolationWindowLowerOffset);
                Console.WriteLine("IsolationWindowUpperOffset: {0}", precursorInfo.IsolationWindowUpperOffset);
                Console.WriteLine("MsLevel: {0}", run.GetMsLevel(scanNum));
            }

            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestXCaliburReader()
        {
            var xcaliburReader = new XCaliburReader(TestTopDownRawFilePathCid);
            var scans = new[] {423, 425};
            foreach (var scan in scans)
            {
                var spec = xcaliburReader.ReadMassSpectrum(scan) as ProductSpectrum;
                Assert.True(spec != null);
                var isolationWindow = spec.IsolationWindow;
                Console.WriteLine("MsLevel: {0}", spec.MsLevel);
                Console.WriteLine("ActivationMethod: {0}", spec.ActivationMethod);
                Console.WriteLine("Rt: {0}", spec.ElutionTime);
                Console.WriteLine("IsolationWindowTargetMz: {0}", isolationWindow.IsolationWindowTargetMz);
                Console.WriteLine("IsolationWindowLowerOffset: {0}", isolationWindow.IsolationWindowLowerOffset);
                Console.WriteLine("IsolationWindowUpperOffset: {0}", isolationWindow.IsolationWindowUpperOffset);
                Console.WriteLine("MonoisotopicMz: {0}", isolationWindow.MonoisotopicMz);
                Console.WriteLine("PrecursorCharge: {0}", isolationWindow.Charge);
            }
        }

        [Test]
        public void TestReadingDiaRawFile()
        {
            const string rawFilePath = @"H:\Research\Jarret\20mz\raw\Q_2014_0523_28_100_amol_uL_20mz.raw";
            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(100);
            spec.Display();
        }

        //[Test]
        //public void TestGeneratingProductXic()
        //{
        //    const string rawFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raw";
        //    var run1 = InMemoryLcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);
        //    var run2 = new DiaLcMsRun(new OldPbfReader(Path.ChangeExtension(rawFilePath, ".pbf")), 0.0, 0.0);

        //    var mz = 1401.643;
        //    var tolerance = new Tolerance(10);
        //    var xic1 = run1.GetFullProductExtractedIonChromatogram(mz, tolerance, 791.03);
        //    xic1.Display();
        //    var xic2 = run2.GetFullProductExtractedIonChromatogram(mz, tolerance, 791.03);
        //    //xic2.Display();
        //    Assert.True(xic1.Count == xic2.Count);
        //    for (var i = 0; i < xic1.Count; i++)
        //    {
        //        if (!xic1[i].Equals(xic2[i]))
        //        {
        //            Console.WriteLine("{0} {1} {2}", i, xic1[i], xic2[i]);
        //        }
        //        Assert.True(xic1[i].Equals(xic2[i]));
        //    }
        //    Console.WriteLine("Done");
        //}

        [Test]
        public void TestGeneratingProductManyXics()
        {
            const string rawFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raw";
            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);
            //var run2 = new DiaLcMsRun(new OldPbfReader(Path.ChangeExtension(rawFilePath, ".pbf")), 0.0, 0.0);

            var tolerance = new Tolerance(10);

            var mzArr = new double[100000];
            var precursorMzArr = new double[mzArr.Length];
            var rnd = new Random();
            for (var i = 0; i < mzArr.Length; i++)
            {
                mzArr[i] = rnd.NextDouble()*1450.0 + 50.0;
                precursorMzArr[i] = rnd.NextDouble()*(810.0-390.0) + 390.0;
            }

            var sw = new System.Diagnostics.Stopwatch();
            double sec;

            // method 1
            sw.Start();
            for (var i = 0; i < mzArr.Length; i++)
            {
                var mz = mzArr[i];
                var tolTh = tolerance.GetToleranceAsTh(mz);
                var minMz = mz - tolTh;
                var maxMz = mz + tolTh;
                var xic1 = run.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
                //var xic2 = run.GetFullProductExtractedIonChromatogram2(minMz, maxMz, precursorMzArr[i]);
                //Assert.True(xic1.Equals(xic2));
            }
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Method 1: {0:f4} sec", sec);

            sw.Reset();
            sw.Start();
            for (var i = 0; i < mzArr.Length; i++)
            {
                var mz = mzArr[i];
                var tolTh = tolerance.GetToleranceAsTh(mz);
                var minMz = mz - tolTh;
                var maxMz = mz + tolTh;
                run.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
            }
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Method 2: {0:f4} sec", sec);

            Console.WriteLine("Done");
        }

        [Test]
        public void TestNoiseFiltration()
        {
            //var run = new PbfLcMsRun(TestQExactiveRawFilePath);
            var run = PbfLcMsRun.GetLcMsRun(TestQExactiveRawFilePath, MassSpecDataType.XCaliburRun);
            //var run = InMemoryLcMsRun.GetLcMsRun(TestQExactiveRawFilePath, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(6747);
            var filtered = spec.GetFilteredSpectrumBySignalToNoiseRatio();
            var productSpec = filtered as ProductSpectrum;
            Console.WriteLine(productSpec == null);
        }

        [Test]
        public void TestReadingCorruptedRawFile()
        {
            const string rawFilePath = @"H:\Research\Corrupted\YS_Shew_testHCD_CID.raw";
            var reader = new XCaliburReader(rawFilePath);
            reader.ReadMassSpectrum(17957);
            try
            {
                reader.ReadMassSpectrum(17957);
            }
            catch (System.Runtime.InteropServices.COMException e)
            {
                Console.WriteLine(e.Message);
            }

            Console.Write("Done");
        }
    }
}
