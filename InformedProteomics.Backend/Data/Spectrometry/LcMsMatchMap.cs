using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class LcMsMatchMap
    {
        public LcMsMatchMap()
        {
            _map = new Dictionary<int, IList<IntRange>>();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass, Tolerance tolerance, LcMsRun run)
        {
            var massBinNum = GetBinNumber(sequenceMass);
            IEnumerable<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(massBinNum, out ms2ScanNums)) return ms2ScanNums;
            return new int[0];
        }

        public void CreateSequenceMassToMs2ScansMap(LcMsRun run, Tolerance tolerance, double minMass, double maxMass)
        {
            // Make a bin to scan numbers map without considering tolerance
            var massBinToScanNumsMapNoTolerance = new Dictionary<int, List<int>>();
            var minBinNum = GetBinNumber(minMass);
            var maxBinNum = GetBinNumber(maxMass);
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                IList<IntRange> scanRanges;
                if (!_map.TryGetValue(binNum, out scanRanges)) continue;
                var sequenceMass = GetMass(binNum);
                var ms2ScanNums = new List<int>();

                foreach (var scanRange in scanRanges)
                {
                    for (var scanNum = scanRange.Min; scanNum <= scanRange.Max; scanNum++)
                    {
                        
                        if (scanNum < run.MinLcScan || scanNum > run.MaxLcScan) continue;
                        if (run.GetMsLevel(scanNum) == 2)
                        {
                            var productSpec = run.GetSpectrum(scanNum) as ProductSpectrum;
                            if (productSpec == null) continue;
                            var isolationWindow = productSpec.IsolationWindow;
                            var isolationWindowTargetMz = isolationWindow.IsolationWindowTargetMz;
                            var charge = (int)Math.Round(sequenceMass / isolationWindowTargetMz);
                            var mz = Ion.GetIsotopeMz(sequenceMass, charge,
                                Averagine.GetIsotopomerEnvelope(sequenceMass).MostAbundantIsotopeIndex);
                            if (productSpec.IsolationWindow.Contains(mz)) ms2ScanNums.Add(scanNum);
                        }
                    }
                }

                ms2ScanNums.Sort();
                massBinToScanNumsMapNoTolerance.Add(binNum, ms2ScanNums);
            }

            // Account for mass tolerance
            _sequenceMassBinToScanNumsMap = new Dictionary<int, IEnumerable<int>>();
            var sumScanNums = 0L;
            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                var sequenceMass = GetMass(binNum);
                var deltaMass = tolerance.GetToleranceAsDa(sequenceMass, 1);

                var curMinBinNum = GetBinNumber(sequenceMass - deltaMass);
                var curMaxBinNum = GetBinNumber(sequenceMass + deltaMass);

                var ms2ScanNums = new HashSet<int>();
                for (var curBinNum = curMinBinNum; curBinNum <= curMaxBinNum; curBinNum++)
                {
                    if (curBinNum < minBinNum || curBinNum > maxBinNum) continue;
                    List<int> existingMs2ScanNums;
                    if (!massBinToScanNumsMapNoTolerance.TryGetValue(curBinNum, out existingMs2ScanNums)) continue;
                    foreach (var ms2ScanNum in existingMs2ScanNums)
                    {
                        ms2ScanNums.Add(ms2ScanNum);
                    }
                }
                _sequenceMassBinToScanNumsMap[binNum] = ms2ScanNums.ToArray();
                sumScanNums += ms2ScanNums.Count;
            }
            Console.WriteLine("#MS/MS matches per sequence: {0}", sumScanNums / (float)(maxBinNum-minBinNum+1));
            _map = null;
        }

        public void SetMatches(double monoIsotopicMass, int minScanNum, int maxScanNum)
        {
            var range = new IntRange(minScanNum, maxScanNum);

            var binNum = GetBinNumber(monoIsotopicMass);
            IList<IntRange> ranges;
            if (_map.TryGetValue(binNum, out ranges))
            {
                var newRanges = new List<IntRange>();
                foreach (var existingRange in ranges)
                {
                    if (range.Overlaps(existingRange))
                    {
                        range = IntRange.Union(range, existingRange);
                    }
                    else
                    {
                        newRanges.Add(existingRange);
                    }
                }
                newRanges.Add(range);
                _map[binNum] = newRanges;
            }
            else
            {
                _map[binNum] = new List<IntRange> {range};
            }
        }

        public static int GetBinNumber(double mass)
        {
            return (int)Math.Round(mass * Constants.RescalingConstantHighPrecision);
        }

        public static double GetMass(int binNum)
        {
            return binNum / Constants.RescalingConstantHighPrecision;
        }

        private Dictionary<int, IList<IntRange>> _map;
        private new Dictionary<int, IEnumerable<int>> _sequenceMassBinToScanNumsMap;
    }
}
