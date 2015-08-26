using System;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Reflection;
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
    public class TestPbfLcMsRun
    {
        const string TestRawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";

        [Test]
        public void TestWritingPbfFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestRawFilePath);
                return;
            }

            Console.WriteLine(@"Writing...");
            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
  
            var outputFilePath = PbfLcMsRun.GetPbfFileName(TestRawFilePath);
            var pbf = new PbfLcMsRun(TestRawFilePath, null, outputFilePath);
            
            Console.WriteLine(@"Done. {0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestReadingPbfFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string pbfFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.pbf";
            if (!File.Exists(pbfFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, pbfFilePath);
                return;
            }

            var pbfRun = new PbfLcMsRun(pbfFilePath);

            var specFilePath = Path.ChangeExtension(pbfFilePath, "raw");
            if (!File.Exists(specFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, specFilePath);
                return;
            }

            Console.WriteLine(@"Loading .pbf into memory");

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath);

            Console.WriteLine(@"Comparing spectra between .pbf and in-memory spectra");

            // spectrum comparison
            for (var scanNum = run.MinLcScan; scanNum <= run.MaxLcScan; scanNum++)
            {
                var spec1 = run.GetSpectrum(scanNum);
                var spec2 = pbfRun.GetSpectrum(scanNum);

                Assert.IsTrue(spec1.Peaks.Length == spec2.Peaks.Length);
                for (var i = 0; i < spec1.Peaks.Length; i++)
                {
                    var p1 = spec1.Peaks[i];
                    var p2 = spec2.Peaks[i];

                    Assert.True(p1.Equals(p2));

                    Assert.True(Math.Abs(p1.Mz - p2.Mz) < 1e-8);
                    Assert.True(Math.Abs(p1.Intensity - p2.Intensity) < 0.001);
                }
            }

            Console.WriteLine(@"Comparing XICs");
            // chromatogram comparison
            const double targetMz = 655.01;
            var tolerance = new Tolerance(10);
            var xic1 = run.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);
            var xic2 = pbfRun.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);
            Assert.True(xic1.Count == xic2.Count);

            for (var i = 0; i < xic1.Count; i++)
            {
                if (!xic1[i].Equals(xic2[i]))
                {
                    Console.WriteLine(@"{0} {1} {2}", i, xic1[i], xic2[i]);
                }
                Assert.True(xic1[i].Equals(xic2[i]));
            }
            Console.WriteLine(@"Done");
        }

        [Test]
        public void TestPbf()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var jobid = new int[]
            {
                1178207, 1178208, 1178209, 1178210, 1178211, 1178212, 1178213, 1178214, 1178215, 1178216, 1178217,
                1178218, 1178219, 1178220, 1178221,
                1178222, 1178223, 1178224, 1178225, 1178226, 1178227, 1178228, 1178229, 1178230, 1178231, 1178232,
                1178233, 1178234, 1178235, 1178236,
                1178237, 1178238, 1178239, 1178240, 1178241, 1178242, 1178243, 1178244, 1178245, 1178246, 1178247,
                1178248, 1178249, 1178250, 1178251,
                1178252, 1178253, 1178254, 1178255, 1178256, 1178257
            };

            var did = new int[]
            {
                02, 01, 04, 05, 03, 07, 10, 08, 09, 11, 14, 17, 25, 21, 24, 15, 26,
                20, 19, 22, 16, 12, 06, 18, 28, 23, 27, 29, 36, 30, 37, 38, 33,
                46, 40, 45, 44, 41, 34, 13, 43, 31, 32, 39, 42, 47, 35, 50, 48, 51, 49
            };
            
            //const string pbfPath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2014_3";
            //const string ms1ftPAth = @"\\proto-5\VOrbiETD02\2014_3\Lewy_intact_02\PMX201503271049_Auto{0}";

            for (var i = 0; i < 51; i++)
            {

                var id = string.Format("00{0}", did[i]);
                id = id.Substring(id.Length - 2);
                var j = jobid[i];

                var srcFile = string.Format(
                    @"\\proto-5\VOrbiETD02\2014_3\Lewy_intact_{0}\PMX201503271049_Auto{1}\Lewy_intact_{2}.ms1ft", id, j, id);

                if (!File.Exists(srcFile))
                {
                    Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, srcFile);
                    return;
                }

                //var destFile = string.Format(@"D:\MassSpecFiles\Lewy\Lewy_intact_{0}.ms1ft", id);
                //File.Copy(srcFile, destFile);
                //var fname = string.Format(@"{0}\Lewy_intact_{1}.pbf", pbfPath, id);
                //var pbfRun = new PbfLcMsRun(fname);
                //Console.WriteLine("{0}\t{1}", pbfRun.MaxLcScan, pbfRun.GetElutionTime(pbfRun.MaxLcScan));
            }
        }

        [Test]
        public void TestGetChrom()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rafFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.pbf";
            if (!File.Exists(rafFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rafFilePath);
                return;
            }

            var pbfRun = new PbfLcMsRun(rafFilePath);

            Console.WriteLine(@"Chromatogram");
            // chromatogram comparison
            const double targetMz = 655.01;
            var tolerance = new Tolerance(10);
            var xic = pbfRun.GetFullPrecursorIonExtractedIonChromatogram(targetMz, tolerance);

            const int MAX_POINTS = 50;
            xic.Display(MAX_POINTS);

            Console.WriteLine(@"Done");
        }

        [Test]
        public void TestRunningTimeChromGen()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            if (!File.Exists(rafFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rafFilePath);
                return;
            }

            var rafRun = new PbfLcMsRun(rafFilePath);

            var tolerance = new Tolerance(10);

            const string dbFile = @"D:\Research\Data\CommonContaminants\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            if (!File.Exists(dbFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, dbFile);
                return;
            }

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
            
            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestGeneratingProductXic()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestRawFilePath);
                return;
            }

            var run = InMemoryLcMsRun.GetLcMsRun(TestRawFilePath);

//            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            if (!File.Exists(rafFilePath))
            {
                Console.WriteLine(@"Warning: Skipping raf portion of test {0} since file not found: {1}", methodName, rafFilePath);
                return;
            }

            var rafRun = new PbfLcMsRun(rafFilePath);

            const double precursorIonMz = 815.16;
            const double productIonMz = 902.445;
            var tolerance = new Tolerance(10);
            var xic1 = run.GetFullProductExtractedIonChromatogram(productIonMz, tolerance, precursorIonMz);
//            xic1.Display();
            var xic2 = rafRun.GetFullProductExtractedIonChromatogram(productIonMz, tolerance, precursorIonMz);
//            xic2.Display();
            Assert.True(xic1.Equals(xic2));
            Console.WriteLine(@"Done");
        }

        [Test]
        public void TestGeneratingProductXics()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestRawFilePath);
                return;
            }

            var run = InMemoryLcMsRun.GetLcMsRun(TestRawFilePath);

//            const string rafFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raf";
            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            if (!File.Exists(rafFilePath))
            {
                Console.WriteLine(@"Warning: Skipping raf portion of test {0} since file not found: {1}", methodName, rafFilePath);
                return;
            }

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
            //double sec;

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
            
            Console.WriteLine(@"Method 1: {0:f4} sec", sw.Elapsed.TotalSeconds);

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
            
            Console.WriteLine(@"Method 2: {0:f4} sec", sw.Elapsed.TotalSeconds);

            Console.WriteLine(@"Done");
        }

        [Test]
        public void TestGeneringAllXics()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(TestRawFilePath, 0.0, 0.0);

            Assert.True(run != null);

            var ms1Scans = run.GetMs1ScanVector();

            var comparer = new MzComparerWithBinning(27);
            const double minMz = 600.0; // 600.0
            const double maxMz = 2000.0;    // 2000.0
            var minBinNum = comparer.GetBinNumber(minMz);
            var maxBinNum = comparer.GetBinNumber(maxMz);
            var scanCount = run.MaxLcScan - run.MinLcScan + 1;

            var sw = new Stopwatch();
            sw.Start();
            var numBinsProcessed = 0;
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
               
                var mzStart = comparer.GetMzStart(binNum);
                var mzEnd = comparer.GetMzEnd(binNum);
                
                var vec1 = run.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);

                if (++numBinsProcessed % 2500 == 0)
                {
                    Console.WriteLine(@"{0,5} / {1,5}, {2:0.00} m/z", numBinsProcessed, maxBinNum - minBinNum + 1, mzStart);
                }
                Assert.IsTrue(vec1.Length == ms1Scans.Length);
            }
            sw.Stop();

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestGettingXicVector()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test " + methodName + @" since file not found: " + TestRawFilePath);
                return;
            }

            var run1 = PbfLcMsRun.GetLcMsRun(TestRawFilePath, 0.0, 0.0);
            var run2 = InMemoryLcMsRun.GetLcMsRun(TestRawFilePath, 0.0, 0.0);

            Assert.True(run1 != null && run2 != null);
            var comparer = new MzComparerWithBinning(27);
            const double minMz = 600.0; // 600.0
            const double maxMz = 2000.0;    // 2000.0
            var minBinNum = comparer.GetBinNumber(minMz);
            var maxBinNum = comparer.GetBinNumber(maxMz);
            Console.WriteLine(@"NumBins: " + (maxBinNum - minBinNum));

            var warnCount = 0;

            var sw = new Stopwatch();
            sw.Start();
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                var mzStart = comparer.GetMzStart(binNum);
                var mzEnd = comparer.GetMzEnd(binNum);

                var vec1 = run1.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);
                var vec2 = run2.GetFullPrecursorIonExtractedIonChromatogramVector(mzStart, mzEnd);

                Assert.True(vec1.Length == vec2.Length,
                                "Extracted Ion Chromatogram vector length mismatch, {0} vs. {1}", vec1.Length, vec2.Length);

                var matchCount = 0;
                var misMatchCount = 0;
                var positiveCount = 0;

                for (var i = 0; i < vec2.Length; i++)
                {
                    if (vec1[i] > 0 || vec2[i] > 0)
                    {
                        positiveCount++;
                    }

                    if (Math.Abs(vec1[i] - vec2[i]) < float.Epsilon)
                    {
                        matchCount++;
                    }
                    else
                    {
                        misMatchCount++;
                    }
                                               
                }

                if (misMatchCount <= 0 || positiveCount == 0)
                {
                    continue;
                }

                var fractionAgreement = 1 - misMatchCount / (double)positiveCount;

                if (fractionAgreement < 0.80 && misMatchCount > 2)
                {
                    Console.WriteLine(@"{0}/{1} Xic values do not match for bin {2} ({3:0.00} m/z); {4:0.0}% agreement",
                                      misMatchCount, positiveCount, binNum, mzStart, fractionAgreement * 100);
                    warnCount++;
                }
            }
            sw.Stop();

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);

            Assert.IsTrue(warnCount < 10, "Too many Xic mismatch warnings: {0}", warnCount);

        }

        [Test]
        public void TestSpectrumNavigation()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rafFilePath = @"H:\Research\Jarret\10mz\raw\Q_2014_0523_50_10_fmol_uL_10mz.raf";
            if (!File.Exists(rafFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rafFilePath);
                return;
            }

            var rafRun = new PbfLcMsRun(rafFilePath);
           // var spec = rafRun.GetSpectrum(54531);
            var scanNum = rafRun.GetNextScanNum(54530, 1);
            Console.WriteLine(scanNum);
        }

        [Test]
        public void TestPbfGen()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            if (!File.Exists(TestRawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, TestRawFilePath);
                return;
            }

            var args = new[] { "-s", TestRawFilePath };
            PbfGen.Program.Main(args);
        }
    }
}
