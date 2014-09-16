using System;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestRafLcMsRun
    {
        [Test]
        public void TestWritingRafFile()
        {
//            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            //const string specFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raw";
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownYufeng\raw\yufeng_column_test2.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            Console.WriteLine("Writing...");
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
  //          const string outputFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
//            const string outputFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            var outputFilePath = Path.ChangeExtension(specFilePath, ".raf");
            run.WriteAsPbf(outputFilePath);
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestReadingRafFile()
        {
            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);

            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            // spectrum comparison
            for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            {
                var spec1 = run.GetSpectrum(scanNum);
                var spec2 = rafRun.GetSpectrum(scanNum);

                Assert.IsTrue(spec1.Peaks.Length == spec2.Peaks.Length);
                for (var i = 0; i < spec1.Peaks.Length; i++)
                {
                    var p1 = spec1.Peaks[i];
                    var p2 = spec2.Peaks[i];
                    Assert.True(p1.Equals(p2));
                }
            }

            Console.WriteLine("Chromatogram");
            // chromatogram comparison
            const double targetMz = 655.01;
            var tolerance = new Tolerance(10);
            var xic1 = run.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);
            var xic2 = rafRun.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);
            Assert.True(xic1.Count == xic2.Count);
            for (var i = 0; i < xic1.Count; i++)
            {
                if (!xic1[i].Equals(xic2[i]))
                {
                    Console.WriteLine("{0} {1} {2}", i, xic1[i], xic2[i]);
                }
                Assert.True(xic1[i].Equals(xic2[i]));
            }
            Console.WriteLine("Done");
        }

        [Test]
        public void TestGetChrom()
        {
            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);

            Console.WriteLine("Chromatogram");
            // chromatogram comparison
            const double targetMz = 655.01;
            var tolerance = new Tolerance(10);
            var xic = rafRun.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);
            xic.Display();

            Console.WriteLine("Done");
        }

        [Test]
        public void TestRunningTimeChromGen()
        {
            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            //var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);

            var tolerance = new Tolerance(10);

            const string dbFile = @"D:\Research\Data\CommonContaminants\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            var numPeptides = 0;
            foreach (var peptide in indexedDb.AnnotationsAndOffsets(6, 30, 2, 2, Enzyme.Trypsin))
            {
                ++numPeptides;
                var comp = new Sequence(peptide.Annotation.Substring(2, peptide.Annotation.Length-4), aaSet).Composition + Composition.H2O;
                var mz = new Ion(comp, 2).GetMonoIsotopicMz();
                //Console.WriteLine(peptide.Annotation + " " + mz);
                rafRun.GetFullPrecursorIonExtractedIonChromatogram(mz, tolerance);
                //run.GetFullPrecursorIonExtractedIonChromatogram(mz, tolerance);

                //var xic1 = run.GetFullPrecursorIonExtractedIonChromatogram(mz, tolerance);
                //var xic2 = rafRun.GetFullPrecursorIonExtractedIonChromatogram(mz, tolerance);
                //Assert.True(xic1.Count == xic2.Count);
                //for (var i = 0; i < xic1.Count; i++)
                //{
                //    if (!xic1[i].Equals(xic2[i]))
                //    {
                //        Console.WriteLine("{0} {1} {2}", i, xic1[i], xic2[i]);
                //    }
                //    Assert.True(xic1[i].Equals(xic2[i]));
                //}

                if (numPeptides == 100000) break;
            }
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestGeneratingProductXic()
        {
//            const string rawFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            const string rawFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raw";

            var run = LcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);

//            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);

            const double precursorIonMz = 815.16;
            const double productIonMz = 902.445;
            var tolerance = new Tolerance(10);
            var xic1 = run.GetFullProductExtractedIonChromatogram(productIonMz, tolerance, precursorIonMz);
//            xic1.Display();
            var xic2 = rafRun.GetFullProductExtractedIonChromatogram(productIonMz, tolerance, precursorIonMz);
//            xic2.Display();
            Assert.True(xic1.Equals(xic2));
            Console.WriteLine("Done");
        }

        [Test]
        public void TestGeneratingProductXics()
        {
//            const string rawFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            const string rawFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raw";
            var run = LcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);

//            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);

            var tolerance = new Tolerance(10);

            var mzArr = new double[100000];
            var precursorMzArr = new double[mzArr.Length];
            var rnd = new Random();
            for (var i = 0; i < mzArr.Length; i++)
            {
                mzArr[i] = rnd.NextDouble() * 1450.0 + 50.0;
                precursorMzArr[i] = rnd.NextDouble() * (810.0 - 390.0) + 390.0;
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
                //var xic2 = rafRun.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
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
                rafRun.GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorMzArr[i]);
            }
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Method 2: {0:f4} sec", sec);

            Console.WriteLine("Done");
        }

        [Test]
        public void TestSpectrumNavigation()
        {
            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            var rafRun = new PbfLcMsRun(rafFilePath);
           // var spec = rafRun.GetSpectrum(54531);
            var scanNum = rafRun.GetNextScanNum(54530, 1);
            Console.WriteLine(scanNum);
        }

        [Test]
        public void TestPbfGen()
        {
            const string rafFilePath = @"H:\Research\Yufeng\TopDownYufeng\raw\yufeng_column_test2.raw";
            var args = new[] {"-s", rafFilePath};
            PbfGen.Program.Main(args);
        }
    }
}
