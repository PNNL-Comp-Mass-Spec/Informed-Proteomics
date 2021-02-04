using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class MsDeconvFilter : ISequenceFilter
    {
        public MsDeconvFilter(LcMsRun run, Tolerance massTolerance, string msDeconvFileName)
        {
            _run = run;
            _massTolerance = massTolerance;
            _lcMsMatchMap = new LcMsMatchMap();

            Read(msDeconvFileName);
        }

        private readonly LcMsRun _run;
        private readonly Tolerance _massTolerance;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass, _massTolerance, _run);
        }

        private void Read(string msDeconvFileName)
        {
            var curScan = 0;
            var isMs1 = false;
            var minMass = double.MaxValue;
            var maxMass = 0.0;

            var featureCountUnfiltered = 0;
            var featureCountFiltered = 0;

            foreach (var line in File.ReadLines(msDeconvFileName))
            {
                if (line.StartsWith("SCANS="))
                {
                    curScan = Convert.ToInt32(line.Substring(line.LastIndexOf('=') + 1));
                    isMs1 = _run.GetMsLevel(curScan) == 1;
                }
                else if(isMs1 && line.Length > 0)
                {
                    var token = line.Split();
                    if (token.Length != 3)
                    {
                        continue;
                    }

                    featureCountUnfiltered++;

                    var charge = Convert.ToInt32(token[2]);
                    if (charge == 1)
                    {
                        continue;
                    }

                    var monoMass = Convert.ToDouble(token[0]);

                    if (minMass > monoMass)
                    {
                        minMass = monoMass;
                    }

                    if (maxMass < monoMass)
                    {
                        maxMass = monoMass;
                    }

                    featureCountFiltered++;

                    var minScan = _run.GetPrevScanNum(curScan, 1);
                    var maxScan = _run.GetNextScanNum(curScan, 1);
                    _lcMsMatchMap.SetMatches(monoMass, minScan, maxScan);
                }
            }

            Console.Write(@"{0}/{1} features loaded...", featureCountFiltered, featureCountUnfiltered);

            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _massTolerance, minMass, maxMass);
        }
    }
}
