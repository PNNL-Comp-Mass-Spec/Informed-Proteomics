using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.UnitTests
{
    [TestFixture]
    public class TestLcMsRun
    {
        public void TestReadingIsolationWindows()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(FilePaths.TestRawFilePath))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since file not found: " + FilePaths.TestRawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(FilePaths.TestRawFilePath, 10000, 10100);
            //for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            foreach (var scanNum in run.AllScanNumbers)
            {
                var isolationWindow = run.GetIsolationWindow(scanNum);
                if (isolationWindow != null)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", scanNum, isolationWindow.MonoisotopicMz ?? 0.0, isolationWindow.Charge ?? 0.0);
                }
            }
        }

        public void TestReadingScanNums()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(FilePaths.TestRawFilePath))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since file not found: " + FilePaths.TestRawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(FilePaths.TestRawFilePath, 20000, 20100);

            var msLevel = new Dictionary<int, int>();

            //for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            foreach (var scanNum in run.AllScanNumbers)
            {
                msLevel[scanNum] = run.GetMsLevel(scanNum);
            }

            //for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            foreach (var scanNum in run.AllScanNumbers)
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
                    //for (var nextScan = scanNum + 1; nextScan <= run.MaxLcScan; nextScan++)
                    foreach (var nextScan in run.AllScanNumbers.Where(x => x > scanNum))
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
            Assert.True(run.GetNextScanNum(20025) == 20032);
            Console.WriteLine(run.GetNextScanNum(20025));
        }

        public void TestParsingSpectrumFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            if (!File.Exists(FilePaths.TestTopDownRawFilePathCid))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, FilePaths.TestTopDownRawFilePathCid);
            }

            const int SCAN = 425;
            const int MAX_POINTS = 50;

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(FilePaths.TestTopDownRawFilePathCid, SCAN, SCAN);

            const int scanNum = SCAN;
            var spec = run.GetSpectrum(scanNum) as ProductSpectrum;
            if (spec != null)
            {
                spec.Display(MAX_POINTS);
                var precursorInfo = spec.IsolationWindow;
                Console.WriteLine("ActivationMethod: {0}", spec.ActivationMethod);
                Console.WriteLine("Rt: {0}", spec.ElutionTime);
                Console.WriteLine("PrecursorScan: {0}", run.GetPrecursorScanNum(spec.ScanNum));
                Console.WriteLine("IsolationWindowTargetMz: {0}", precursorInfo.IsolationWindowTargetMz);
                Console.WriteLine("IsolationWindowLowerOffset: {0}", precursorInfo.IsolationWindowLowerOffset);
                Console.WriteLine("IsolationWindowUpperOffset: {0}", precursorInfo.IsolationWindowUpperOffset);
                Console.WriteLine("MsLevel: {0}", run.GetMsLevel(scanNum));
            }

            Console.WriteLine(@"Done. {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestXCaliburReader()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var xcaliburReader = new XCaliburReader(FilePaths.TestTopDownRawFilePathCid);
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
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = FilePaths.TestRawFilePath;
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            const int SCAN = 100;
            const int MAX_POINTS = 50;

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(rawFilePath, SCAN);
            var spec = run.GetSpectrum(SCAN);
            spec.Display(MAX_POINTS);

            // Console.WriteLine("{0}, {1}", spec.Peaks[50].Mz, spec.Peaks[50].Intensity);
            // Console.WriteLine("{0}, {1}", spec.Peaks[500].Mz, spec.Peaks[500].Intensity);
            // Console.WriteLine("{0}, {1}", spec.Peaks[1000].Mz, spec.Peaks[1000].Intensity);

            Assert.IsTrue(Math.Abs(spec.Peaks[50].Mz - 414.75503540039062) < 0.0001, "Invalid m/z for peak at index 50");
            Assert.IsTrue(Math.Abs(spec.Peaks[50].Intensity - 1071.5673828125) < 0.01, "Invalid intensity for peak at index 50");

            Assert.IsTrue(Math.Abs(spec.Peaks[500].Mz - 578.1298828125) < 0.0001, "Invalid m/z for peak at index 500");
            Assert.IsTrue(Math.Abs(spec.Peaks[500].Intensity - 573.02374267578125) < 0.01, "Invalid intensity for peak at index 500");

            Assert.IsTrue(Math.Abs(spec.Peaks[1000].Mz - 974.17694091796875) < 0.0001, "Invalid m/z for peak at index 1000");
            Assert.IsTrue(Math.Abs(spec.Peaks[1000].Intensity - 678.13824462890625) < 0.01, "Invalid intensity for peak at index 1000");
        }

        //[Test]
        //public void TestGeneratingProductXic()
        //{
        //    var methodName = MethodBase.GetCurrentMethod().Name;
        //    TestUtils.ShowStarting(methodName);

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
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = FilePaths.TestRawFilePath;
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath);
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
            //double sec;

            // method 1
            sw.Start();
            for (var i = 0; i < mzArr.Length; i++)
            {
                var mz = mzArr[i];
                var tolTh = tolerance.GetToleranceAsMz(mz);
                var minMz = mz - tolTh;
                var maxMz = mz + tolTh;
                var xic1 = run.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
                //var xic2 = run.GetFullProductExtractedIonChromatogram2(minMz, maxMz, precursorMzArr[i]);
                //Assert.True(xic1.Equals(xic2));
            }
            sw.Stop();

            Console.WriteLine(@"Method 1: {0:f4} sec", sw.Elapsed.TotalSeconds);

            sw.Reset();
            sw.Start();
            for (var i = 0; i < mzArr.Length; i++)
            {
                var mz = mzArr[i];
                var tolTh = tolerance.GetToleranceAsMz(mz);
                var minMz = mz - tolTh;
                var maxMz = mz + tolTh;
                run.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
            }
            sw.Stop();

            Console.WriteLine(@"Method 2: {0:f4} sec", sw.Elapsed.TotalSeconds);

            Console.WriteLine("Done");
        }

        [Test]
        public void TestNoiseFiltration()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(FilePaths.TestQExactiveRawFilePath))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since file not found: " + FilePaths.TestQExactiveRawFilePath);
            }

            //var run = new PbfLcMsRun(TestQExactiveRawFilePath);
            var run = PbfLcMsRun.GetLcMsRun(FilePaths.TestQExactiveRawFilePath);
            //var run = InMemoryLcMsRun.GetLcMsRun(TestQExactiveRawFilePath, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(6747);
            var filtered = spec.GetFilteredSpectrumBySignalToNoiseRatio();
            var productSpec = filtered as ProductSpectrum;
            Console.WriteLine(productSpec == null);
        }

        [Test]
        [Ignore("Ignore this since it crashes Nunit 3")]
        public void TestReadingCorruptedRawFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\Corrupted\YS_Shew_testHCD_CID.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

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

        [Test]
        public void TestReadingSingleSpecMzMlFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string filePath = @"D:\Research\Data\TRex\VNVADCGAEALAR.mzML";
            if (!File.Exists(filePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, filePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(filePath);
            Console.WriteLine(run.MaxLcScan);
        }

        [Test]
        public void TestReadingBrukerDaltonDataSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"D:\MassSpecFiles\ICR\20141212FGWT1_F5_1_01_3230.mzML";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            Console.WriteLine(@"Test start");

            var pbf = new PbfLcMsRun(specFilePath, null, null, 1.4826, 1.4826);
            //var pbf = new PbfLcMsRun(specFilePath, null, null, 3, 1.4826);
            /*
            foreach (var spec in reader.ReadAllSpectra())
            {
                if (spec != null)
                {
                    Console.WriteLine(spec.ScanNum);
                }
            }*/
            Console.WriteLine(@"Test end");
        }

        [Test]
        public void TestReadingRawFileWithSingleMs2Spectrum()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SpecFiles\2015-05-06_Carbonic_HCD_854_50AVG.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

//            var run = (InMemoryLcMsRun) InMemoryLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);
//            run.WriteAsPbf(Path.ChangeExtension(specFilePath, ".pbf"));
            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
        }
    }
}
