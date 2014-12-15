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
        public Ms1FtFilter(LcMsRun run, Tolerance massTolerance, string isosFileName)
        {
            _run = run;
            _massTolerance = massTolerance;
            _lcMsMatchMap = new LcMsMatchMap();

            Read(isosFileName);
        }

        private readonly LcMsRun _run;
        private readonly Tolerance _massTolerance;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass, _massTolerance, _run);
        }

        private void Read(string ms1FtFileName)
        {
            var ftFileParser = new TsvFileParser(ms1FtFileName);
            var monoMassArr = ftFileParser.GetData("monoisotopic_mw").Select(Convert.ToDouble).ToArray();
            //var repScanArray = ftFileParser.GetData("rep_scan_num").Select(s => Convert.ToInt32(s)).ToArray();
            var minScanArray = ftFileParser.GetData("min_scan_num").Select(s => Convert.ToInt32(s)).ToArray();
            var maxScanArray = ftFileParser.GetData("max_scan_num").Select(s => Convert.ToInt32(s)).ToArray();

            var minMass = double.MaxValue;
            var maxMass = 0.0;
            for (var i = 0; i < monoMassArr.Length; i++)
            {
                //var repScan = repScanArray[i];
                var monoMass = monoMassArr[i];
                if (minMass > monoMass) minMass = monoMass;
                if (maxMass < monoMass) maxMass = monoMass;

                //var minScan = _run.GetPrevScanNum(repScan, 1);
                //var maxScan = _run.GetNextScanNum(repScan, 1);
                var minScan = minScanArray[i];
                var maxScan = maxScanArray[i];
                _lcMsMatchMap.SetMatches(monoMass, minScan, maxScan);
            }

            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _massTolerance, minMass, maxMass);
        }
    }
}
