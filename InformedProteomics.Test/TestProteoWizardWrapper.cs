using System;
using System.Diagnostics;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestProteoWizardWrapper
    {
        public const string TestRawFilePath = @"\\protoapps\UserData\Sangtae\TestData\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";

        [Test]
        public void TestLoadingProteoWizardWrapper()
        {
            //try
            //{
            var reader = new ProteoWizardReader(TestRawFilePath);
            reader.ReadMassSpectrum(1);
            //}
            //catch (FileNotFoundException e)
            //{
            //    Console.WriteLine("FileNotFound: {0}", e.FileName);
            //}
            //var mzs = new[] {1.0};
            //var intensities = new[] {10000.0};
            //var centroider = new Centroider(mzs, intensities);
            //centroider.GetCentroidedData(out mzs, out intensities);
        }

        [Test]
        public void TestReadingXCaliburRun()
        {
            var sw = new Stopwatch();
            sw.Start();
            var run = new LcMsRun(new XCaliburReader(TestRawFilePath));
            Console.WriteLine("{0}\t{1}", run.MinLcScan, run.MaxLcScan);
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }
    }
}
