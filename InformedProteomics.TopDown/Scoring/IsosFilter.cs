using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class IsosFilter : ISequenceFilter
    {
        public IsosFilter(LcMsRun run, Tolerance massTolerance, string isosFileName, double fitScoreThreshold = 1.0)
        {
            _run = run;
            _massTolerance = massTolerance;
            _fitScoreThreshold = fitScoreThreshold;
            _lcMsMatchMap = new LcMsMatchMap();

            Read(isosFileName);
        }

        private readonly LcMsRun _run;
        private readonly Tolerance _massTolerance;
        private readonly double _fitScoreThreshold;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass);
        }

        private void Read(string isosFileName)
        {
            var icrToolsparser = new TsvFileParser(isosFileName, ',');
            var monoMassArr = icrToolsparser.GetData("monoisotopic_mw").Select(Convert.ToDouble).ToArray();
            var scanArray = icrToolsparser.GetData("scan_num").Select(s => Convert.ToInt32(s)).ToArray();
            var chargeArray = icrToolsparser.GetData("charge").Select(s => Convert.ToInt32(s)).ToArray();

            var fitStringArr = icrToolsparser.GetData("fit");
            var fitArray = fitStringArr == null ? null : icrToolsparser.GetData("fit").Select(Convert.ToDouble).ToArray();

            var featureCountFiltered = 0;

            var minMass = double.MaxValue;
            var maxMass = 0.0;
            for (var i = 0; i < monoMassArr.Length; i++)
            {
                if (fitArray != null && fitArray[i] > _fitScoreThreshold || chargeArray[i] <= 1)
                {
                    continue;
                }

                featureCountFiltered++;

                var scan = scanArray[i];
                var monoMass = monoMassArr[i];
                if (minMass > monoMass)
                {
                    minMass = monoMass;
                }

                if (maxMass < monoMass)
                {
                    maxMass = monoMass;
                }

                var minScan = _run.GetPrevScanNum(scan, 1);
                var maxScan = _run.GetNextScanNum(scan, 1);
                _lcMsMatchMap.SetMatches(monoMass, minScan, maxScan);
            }

            Console.Write("{0}/{1} features loaded...", featureCountFiltered, monoMassArr.Length);

            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _massTolerance, minMass, maxMass);
        }
    }
}
