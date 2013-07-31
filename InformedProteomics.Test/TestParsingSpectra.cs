using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.PeakDetectors;
using DeconTools.Backend.Runs;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;
using XCaliburRun = InformedProteomics.Backend.MassSpecData.XCaliburRun;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestParsingSpectra
    {
        public void TestParsingSpectrumFile()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            //var peakDetector = new DeconToolsPeakDetectorV2(0, 0);  // basic peak detector
            //var run = new RunFactory().CreateRun(specFilePath);
            var run = new XCaliburRun(specFilePath);

            //for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            //{
            //    run.GetMassSpectrum(scanNum);
            //}

            const int scanNum = 810;
            var spec = run.GetMassSpectrum(scanNum) as ProductSpectrum;
            if (spec != null)
            {
                spec.Display();
                var precursorInfo = spec.PrecursorInfo;
                Console.WriteLine("ActivationMethod: {0}", spec.ActivationMethod);
                Console.WriteLine("PrecursorScan: {0}", precursorInfo.PrecursorScan);
                Console.WriteLine("IsolationWindowTargetMz: {0}", precursorInfo.IsolationWindowTargetMz);
                Console.WriteLine("IsolationWindowLowerOffset: {0}", precursorInfo.IsolationWindowLowerOffset);
                Console.WriteLine("IsolationWindowUpperOffset: {0}", precursorInfo.IsolationWindowUpperOffset);
                Console.WriteLine("MsLevel: {0}", run.GetMsLevel(scanNum));
                Console.WriteLine("IsCentroidScan: {0}", run.IsCentroidScan(scanNum));
            }

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        public void TestParsingSpectrumFilesUsingDeconTools()
        {
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var peakDetector = new DeconToolsPeakDetectorV2(0, 0);  // basic peak detector
            var run = new RunFactory().CreateRun(specFilePath);

            var numMs2Specs = 0;
            for (var scanNum = run.GetMinPossibleLCScanNum(); scanNum <= run.GetMaxPossibleLCScanNum(); scanNum++)
            {
                if (run.GetMSLevel(scanNum) == 2)
                {
                    ++numMs2Specs;
                    var xyData = run.GetMassSpectrum(new ScanSet(810));
                    peakDetector.FindPeaks(xyData);
                }
            }
            Console.WriteLine("Num MS2 spectra: {0}", numMs2Specs);

            //const int scanNum = 810;
            //var xyData = run.GetMassSpectrum(new ScanSet(scanNum));
            //xyData.Display();
            //peakDetector.FindPeaks(xyData);

            //Console.WriteLine("IsolationWindow: {0}", run.GetMS2IsolationWidth(scanNum));
            //Console.WriteLine("ScanInfo\n{0}", run.GetScanInfo(scanNum));
            //Console.WriteLine("PrecursorScan\n{0}", run.GetPrecursorInfo(scanNum).PrecursorScan); 
            
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }
    }
}
