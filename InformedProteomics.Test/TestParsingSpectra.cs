using System;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.PeakDetectors;
using DeconTools.Backend.Runs;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestParsingSpectra
    {
        [Test]
        public void TestPeakFiltering()
        {
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            run.GetSpectrum(3385).Display();
        }

        [Test]
        public void TestSpecFiltering()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = new LcMsRun(new XCaliburReader(rawFilePath));

            var spec = run.GetSpectrum(3997);
            spec.Display();

            //var numSpecs = 0;
            //var reader = new PgfReader(pgfFilePath);
            //foreach (var spec in reader.ReadAllSpectra())
            //{
            //    ++numSpecs;
            //}
            //Console.WriteLine(@"NumSpectra: {0}", numSpecs);

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestReadingPbf()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            const string pbfFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA_NoCharge.pbf";
            var run = new LcMsRun(new PbfReader(pbfFilePath));
            Console.WriteLine(run.MaxLcScan);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);

            sw.Reset();
            sw.Start();
            const string pgfFilePath = @"D:\Research\Data\UW\QExactive\pgf\82593_lv_mcx_DDA_NoCharge.pgf";
            var run2 = new LcMsRun(new PgfReader(pgfFilePath));
            Console.WriteLine(run.MaxLcScan);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);

        }

        [Test]
        public void TestWritingPgf()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            const string rawFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA.raw";
            var run = new LcMsRun(new XCaliburReader(rawFilePath));

            const string pgfFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA.pgf";
            run.WriteTo(pgfFilePath);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestParsingSpectrumFile()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = new LcMsRun(new XCaliburReader(specFilePath));

            //for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            //{
            //    run.ReadMassSpectrum(scanNum);
            //}

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

        //public void TestParsingSpectrumFilesUsingDeconTools()
        //{
        //    var sw = new System.Diagnostics.Stopwatch();

        //    sw.Start();
            
        //    const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
        //    var peakDetector = new DeconToolsPeakDetectorV2(0, 0);  // basic peak detector
        //    var run = new RunFactory().CreateRun(specFilePath);

        //    var numMs2Specs = 0;
        //    for (var scanNum = run.GetMinPossibleLCScanNum(); scanNum <= run.GetMaxPossibleLCScanNum(); scanNum++)
        //    {
        //        if (run.GetMSLevel(scanNum) == 2)
        //        {
        //            ++numMs2Specs;
        //            var xyData = run.GetMassSpectrum(new ScanSet(810));
        //            peakDetector.FindPeaks(xyData);
        //        }
        //    }
        //    Console.WriteLine("Num MS2 spectra: {0}", numMs2Specs);

        //    //const int scanNum = 810;
        //    //var xyData = run.ReadMassSpectrum(new ScanSet(scanNum));
        //    //xyData.Display();
        //    //peakDetector.FindPeaks(xyData);

        //    //Console.WriteLine("IsolationWindow: {0}", run.GetMS2IsolationWidth(scanNum));
        //    //Console.WriteLine("ScanInfo\n{0}", run.GetScanInfo(scanNum));
        //    //Console.WriteLine("PrecursorScan\n{0}", run.ReadPrecursorInfo(scanNum).PrecursorScan); 
            
        //    var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
        //    Console.WriteLine(@"Done. {0:f4} sec", sec);
        //}
    }
}
