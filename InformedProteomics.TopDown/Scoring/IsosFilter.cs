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
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass, _massTolerance, _run);
        }

        private void Read(string isosFileName)
        {
            var icrToolsparser = new TsvFileParser(isosFileName, ',');
            var monoMassArr = icrToolsparser.GetData("monoisotopic_mw").Select(Convert.ToDouble).ToArray();
            var scanArray = icrToolsparser.GetData("scan_num").Select(s => Convert.ToInt32(s)).ToArray();
            var chargeArray = icrToolsparser.GetData("charge").Select(s => Convert.ToInt32(s)).ToArray();
            var fitArray = icrToolsparser.GetData("fit").Select(Convert.ToDouble).ToArray();

            var minMass = double.MaxValue;
            var maxMass = 0.0;
            for (var i = 0; i < fitArray.Length; i++)
            {
                if (fitArray[i] > _fitScoreThreshold || chargeArray[i] <= 1) continue;
                var scan = scanArray[i];
                var monoMass = monoMassArr[i];
                if (minMass > monoMass) minMass = monoMass;
                if (maxMass < monoMass) maxMass = monoMass;

                var minScan = _run.GetPrevScanNum(scan, 1);
                var maxScan = _run.GetNextScanNum(scan, 1);
                _lcMsMatchMap.SetMatches(monoMass, minScan, maxScan);
            }

            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _massTolerance, minMass, maxMass);
        }
    }
}
