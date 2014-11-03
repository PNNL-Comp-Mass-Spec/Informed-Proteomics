//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using InformedProteomics.Backend.Data.Spectrometry;
//using InformedProteomics.Backend.MassSpecData;
//using InformedProteomics.Backend.Utils;

//namespace InformedProteomics.TopDown.Scoring
//{
//    public class IsosFilter : ISequenceFilter
//    {
//        public IsosFilter(LcMsRun run, Tolerance massTolerance, string isosFileName, double fitScoreThreshold = 1.0)
//        {
//            _run = run;
//            _massTolerance = massTolerance;
//            _fitScoreThreshold = fitScoreThreshold;
//            _lcMsMatchMap = new LcMsMatchMap();

//            Read(isosFileName);
//        }

//        private readonly LcMsRun _run;
//        private readonly Tolerance _massTolerance;
//        private readonly double _fitScoreThreshold;
//        private readonly LcMsMatchMap _lcMsMatchMap;

//        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
//        {
//            var xic = new Xic();

//            var index = peakList.BinarySearch(new LcMsPeak((minMz + maxMz) / 2, 0, 0));
//            if (index < 0) index = ~index;

//            // go down
//            var i = index - 1;
//            while (i >= 0 && i < peakList.Count)
//            {
//                var peak = peakList[i];
//                if (peak.Mz <= minMz) break;
//                xic.Add(new XicPoint(peak.ScanNum, peak.Mz, peak.Intensity));
//                --i;
//            }

//            // go up
//            i = index;
//            while (i >= 0 && i < peakList.Count)
//            {
//                var peak = peakList[i];
//                if (peak.Mz >= maxMz) break;
//                xic.Add(new XicPoint(peak.ScanNum, peak.Mz, peak.Intensity));
//                ++i;
//            }

//            return xic;
//        }

//        private List<Feature> _features;

//        private class Feature
//        {

//        }
//        private void Read(string isosFileName)
//        {
//            var icrToolsparser = new TsvFileParser(isosFileName, ',');
//            var monoMassArr = icrToolsparser.GetData("monoisotopic_mw").Select(Convert.ToDouble).ToArray();
//            var scanArray = icrToolsparser.GetData("scan_num").Select(s => Convert.ToInt32(s)).ToArray();
//            var fitArray = icrToolsparser.GetData("fit").Select(Convert.ToDouble).ToArray();
//            var peakList = fitArray.Where(fit => fit < _fitScoreThreshold).Select((fit, i) => new Peak(monoMassArr[i], scanArray[i])).OrderBy(p => p.Mz).ToList();
//            Console.WriteLine("FitScoreThreshold: {0}", _fitScoreThreshold);
//            Console.WriteLine("NumFeatures: {0}", peakList.Count);
//            //var peakList = monoMassArr.Select((mass, i) => new Peak(mass, scanArray[i])).OrderBy(p => p.Mz).ToList();
//            const double scanTolerance = 10.0;

//            const string resultFilePath = @"H:\Research\Yufeng\TopDownYufeng\M1_V31\yufeng_column_test2_IcTda.tsv";
//            var parser = new TsvFileParser(resultFilePath);
//            var qValues = parser.GetData("QValue").Select(Convert.ToDouble).ToArray();
//            var masses = parser.GetData("Mass").Select(Convert.ToDouble).ToArray();
//            var scans = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
//            var compositions = parser.GetData("Composition").ToArray();
//            var numTotal = 0;
//            var numMatch = 0;
//            var compSet = new HashSet<string>();
//            var matchedCompSet = new HashSet<string>();
//            for (var i = 0; i < qValues.Length; i++)
//            {
//                if (qValues[i] > 0.01) break;
//                var mass = masses[i];
//                var scan = scans[i];
//                var matches = PeakListUtils.FindAllPeaks(peakList, mass, _massTolerance)
//                    .Where(p => p.Intensity >= scan - scanTolerance && p.Intensity <= scan + scanTolerance);
//                if (matches.Any())
//                {
//                    ++numMatch;
//                    matchedCompSet.Add(compositions[i]);
//                }
//                ++numTotal;
//                compSet.Add(compositions[i]);
//            }
//            Console.WriteLine("PrSMs: {0} / {1} = {2}", numMatch, numTotal, numMatch / (float)numTotal);
//            Console.WriteLine("Compositions: {0} / {1} = {2}", matchedCompSet.Count, compSet.Count, matchedCompSet.Count / (float)compSet.Count);
//        }
//    }
//}
