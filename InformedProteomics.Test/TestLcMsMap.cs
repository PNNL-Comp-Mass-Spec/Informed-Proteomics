using System;
using System.Diagnostics;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestLcMsMap
    {
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
