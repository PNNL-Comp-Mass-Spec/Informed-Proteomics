using System;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests
{
    [TestFixture]
    public class TestLcMsMap
    {
        [Test]
        [Category("Local_Testing")]
        public void TestMs1Filter()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // QC_Shew
            const string specFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string ms1FtFileName = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            const string idFilePath =  @"D:\MassSpecFiles\training\IcTda\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            if (!File.Exists(idFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, idFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
            var massTolerance = new Tolerance(10);

            var ms1ftFilter = new Ms1FtFilter(run, massTolerance, ms1FtFileName);
            var n = 0;
            var ms2ScanNums = run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                var matchingMass = ms1ftFilter.GetMatchingMass(ms2ScanNum);
                n += matchingMass.Count();
            }

            Console.WriteLine("{0} / {1}", n, ms2ScanNums.Count);
            /*
            var tsvReader = new TsvFileParser(idFilePath);

            for (var i = 0; i < tsvReader.NumData; i++)
            {
                var qv = double.Parse(tsvReader.GetData("QValue")[i]);
                if (qv > 0.01) break;

                var scan = int.Parse(tsvReader.GetData("Scan")[i]);
                var charge = int.Parse(tsvReader.GetData("Charge")[i]);
                var mass = double.Parse(tsvReader.GetData("Mass")[i]);

                if (mass > 15000) continue;

                var seq = tsvReader.GetData("Sequence")[i];
                var mod = tsvReader.GetData("Modifications")[i];
                var nMatched = int.Parse(tsvReader.GetData("#MatchedFragments")[i]);

                var hit = false;
                foreach (var ms2Scan in ms1ftFilter.GetMatchingMs2ScanNums(mass))
                {
                    if (ms2Scan == scan)
                    {
                        hit = true;
                        break;
                    }
                }

                if (!hit)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", scan, mass, nMatched);
                }
            }*/
        }

        /*
        public void TestCalculatingNumBins()
        {
            const int numBits = 26;
            var mzComparer = new MzComparerWithBinning(numBits);
            const double minMass = 400.0;
            const double maxMass = 4000.0;

            var minBinNum = mzComparer.GetBinNumber(minMass);
            var maxBinNum = mzComparer.GetBinNumber(maxMass);
            var numBins = maxBinNum - minBinNum + 1;

            Console.WriteLine("NumBins: {0}", numBins);
        }

        public void TestRunningTime()
        {
            const string testRawFilePath = @"\\protoapps\UserData\Sangtae\Yufeng\raw\yufeng_column_test2.raw";
            const double monoIsotopicMass = 12377.46179;
            var run = PbfLcMsRun.GetLcMsRun(testRawFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var map = new Ms1FeatureMatrix(run);
            var regions = map.GetProbableChargeScanRegions(monoIsotopicMass);
            var scanNums = map.GetMatchingMs2ScanNums(12377.46179).ToList();
            Console.WriteLine("NumRegions: " + regions.Count());
            Console.WriteLine("NumScanNums: " + scanNums.Count);
            Console.WriteLine(string.Join("\t", scanNums));
        }

        public void TestSummingSpectra()
        {
            //const string testRawFilePath = @"\\protoapps\UserData\Sangtae\Yufeng\raw\yufeng_column_test2.raw";
            const string testRawFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";

            //const double monoIsotopicMass = 12377.46179;

            const int numBits = 26;
            var mzComparer = new MzComparerWithBinning(numBits);
            const double minMass = 3000.0;
            const double maxMass = 50000.0;

            var minBinNum = mzComparer.GetBinNumber(minMass);
            var maxBinNum = mzComparer.GetBinNumber(maxMass);
            var numBins = maxBinNum - minBinNum + 1;
            var sumRegions = 0L;

            Console.WriteLine("NumBins: {0}", numBins);
            Console.Write("Reading raw file...");
            var run = InMemoryLcMsRun.GetLcMsRun(testRawFilePath, MassSpecDataType.XCaliburRun, 0, 1.4826);
            Console.WriteLine(@"Done");

            var map = new Ms1FeatureMatrix(run, numBits);
            var sw = new Stopwatch();
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                if ((binNum - minBinNum)%100 == 0)
                {
                    sw.Stop();
                    Console.WriteLine("{0}% done. {1}", 100*(binNum - minBinNum) / (float)numBins, sw.ElapsedMilliseconds / (float)1000);
                    Console.WriteLine("NumClusters: {0}", sumRegions);
                    sw.Reset();
                    sw.Start();
                }
                var monoIsotopicMass = mzComparer.GetMzAverage(binNum);
                var regions = map.GetProbableChargeScanRegions(monoIsotopicMass);
                foreach (var region in regions)
                {
                    run.GetSummedMs2Spectrum(monoIsotopicMass, region.MinScanNum, region.MaxScanNum, region.MinCharge, region.MaxCharge);
                    ++sumRegions;
                }
            }
            Console.WriteLine("AverageNumRegions: {0}", sumRegions/(float)numBins);
        }

        public void TestOldMs1Filter()
        {
            const string testRawFilePath = @"\\protoapps\UserData\Sangtae\Yufeng\raw\yufeng_column_test2.raw";
            var run = InMemoryLcMsRun.GetLcMsRun(testRawFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var ms1Filter = new Ms1IsotopeAndChargeCorrFilter(run, new Tolerance(15));
            //ms1Filter.PrecomputePossibleSequenceMasses();
            var scanNums = ms1Filter.GetMatchingMs2ScanNums(12377.46179).ToList();
            Console.WriteLine("NumScanNums: " + scanNums.Count);
            Console.WriteLine(string.Join("\t", scanNums));
        }
         */
    }
}
