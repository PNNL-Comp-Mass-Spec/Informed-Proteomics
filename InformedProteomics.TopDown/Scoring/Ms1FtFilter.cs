using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1FtFilter : ISequenceFilter
    {
        public Ms1FtFilter(LcMsRun run, Tolerance massTolerance, string ms1FtFileName, double minLikelihoodRatio = 0)
        {
            Ms1FtIndexToScanRange = new Dictionary<int, Tuple<int, int>>();
            _lcMsChargeMap = new LcMsChargeMap(run, massTolerance);
            _minLikelihoodRatio = minLikelihoodRatio;

            Read(ms1FtFileName);
        }

        private readonly LcMsChargeMap _lcMsChargeMap;

        public Dictionary<int, Tuple<int, int>> Ms1FtIndexToScanRange { get; private set; }

        public IEnumerable<int> GetMatchingFeatureIds(double sequenceMass)
        {
            return _lcMsChargeMap.GetMatchingFeatureIds(sequenceMass);
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsChargeMap.GetMatchingMs2ScanNums(sequenceMass);
        }

        public IEnumerable<double> GetMatchingMass(int ms2Scan)
        {
            return _lcMsChargeMap.GetMatchingMass(ms2Scan);
        }

        private void Read(string ms1FtFileName)
        {
            var ftFileParser = new TsvFileParser(ms1FtFileName);
            var featureIdArr = ftFileParser.GetData("FeatureID").Select(s => Convert.ToInt32(s)).ToArray();
            var monoMassArr = ftFileParser.GetData("MonoMass").Select(Convert.ToDouble).ToArray();
            var minScanArray = ftFileParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            var maxScanArray = ftFileParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            var repScanArray = ftFileParser.GetData("RepScan").Select(s => Convert.ToInt32(s)).ToArray();
            var minChargeArray = ftFileParser.GetData("MinCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var maxChargeArray = ftFileParser.GetData("MaxCharge").Select(s => Convert.ToInt32(s)).ToArray();

            var likelihoodRatioData = ftFileParser.GetData("LikelihoodRatio");

            double[] scoreArray;
            if (likelihoodRatioData != null && likelihoodRatioData.Count > 0)
            {
                scoreArray = likelihoodRatioData.Select(Convert.ToDouble).ToArray();
            }
            else
            {
                scoreArray = new double[monoMassArr.Length];
            }

            var featureCountFiltered = 0;

            for (var i = 0; i < monoMassArr.Length; i++)
            {
                //if (flagArray[i] == 0 && probArray[i] < _minProbability)  continue;
                if (scoreArray[i] < _minLikelihoodRatio) continue;
                featureCountFiltered++;
                var monoMass = monoMassArr[i];
                _lcMsChargeMap.SetMatches(featureIdArr[i], monoMass, minScanArray[i], maxScanArray[i], repScanArray[i], minChargeArray[i], maxChargeArray[i]);
                Ms1FtIndexToScanRange.Add(featureIdArr[i], new Tuple<int, int>(minScanArray[i], maxScanArray[i]));
            }

            // NOTE: The DMS Analysis Manager looks for this statistic; do not change it
            Console.Write(@"{0}/{1} features loaded...", featureCountFiltered, monoMassArr.Length);
            _lcMsChargeMap.CreateMassToScanNumMap();
        }

        private readonly double _minLikelihoodRatio;
    }
}
