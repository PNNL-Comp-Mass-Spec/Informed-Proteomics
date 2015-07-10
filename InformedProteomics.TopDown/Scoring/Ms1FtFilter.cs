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
        public Ms1FtFilter(LcMsRun run, Tolerance massTolerance, string ms1FtFileName)
        {
            _lcMsChargeMap = new LcMsChargeMap(run, massTolerance);

            Read(ms1FtFileName);
        }

        private readonly LcMsChargeMap _lcMsChargeMap;

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsChargeMap.GetMatchingMs2ScanNums(sequenceMass);
        }

        private void Read(string ms1FtFileName)
        {
            var ftFileParser = new TsvFileParser(ms1FtFileName);
            //var featureIdArray = ftFileParser.GetData("FeatureID").Select(s => Convert.ToInt32(s)).ToArray();
            var monoMassArr = ftFileParser.GetData("MonoMass").Select(Convert.ToDouble).ToArray();
            var minScanArray = ftFileParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            var maxScanArray = ftFileParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            var repScanArray = ftFileParser.GetData("RepScan").Select(s => Convert.ToInt32(s)).ToArray();
            var minChargeArray = ftFileParser.GetData("MinCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var maxChargeArray = ftFileParser.GetData("MaxCharge").Select(s => Convert.ToInt32(s)).ToArray();
            //var probArray = ftFileParser.GetData("Probability").Select(Convert.ToDouble).ToArray();
            //var flagArray = ftFileParser.GetData("GoodEnough").Select(s => Convert.ToInt32(s)).ToArray();
            var featureCountFiltered = 0;

            for (var i = 0; i < monoMassArr.Length; i++)
            {
                //if (flagArray[i] == 0 && probArray[i] < _minProbability)  continue;
                featureCountFiltered++;
                var monoMass = monoMassArr[i];
                _lcMsChargeMap.SetMatches(monoMass, minScanArray[i], maxScanArray[i], repScanArray[i], minChargeArray[i], maxChargeArray[i]);
            }

            // NOTE: The DMS Analysis Manager looks for this statistic; do not change it
            Console.Write(@"{0}/{1} features loaded...", featureCountFiltered, monoMassArr.Length);
            _lcMsChargeMap.CreateMassToScanNumMap();
        }
    }
}
