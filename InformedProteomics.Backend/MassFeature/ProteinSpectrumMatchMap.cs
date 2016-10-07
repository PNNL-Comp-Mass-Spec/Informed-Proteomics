using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Integration;

namespace InformedProteomics.Backend.MassFeature
{
    public class ProteinSpectrumMathMap
    {
        public ProteinSpectrumMathMap(LcMsRun run, int dataid, List<ProteinSpectrumMatch> prsmList, string dataDesc = "")
        {
            DataId = dataid;
            DataDesc = dataDesc;
            Run = run;
            _scanNumToMatchMap = new Dictionary<int, ProteinSpectrumMatch>();

            if (prsmList == null) return;

            ProteinSpectrumMatches = prsmList;
            _scanNumToMatchMap.Clear();

            foreach (var prsm in prsmList)
            {
                if (_scanNumToMatchMap.ContainsKey(prsm.ScanNum))
                {
                    if (_scanNumToMatchMap[prsm.ScanNum].Score < prsm.Score) _scanNumToMatchMap[prsm.ScanNum] = prsm;
                }
                else
                {
                    _scanNumToMatchMap.Add(prsm.ScanNum, prsm);
                }
            }
        }

        public ProteinSpectrumMatch FindByFeature(LcMsFeature feature, Tolerance tolerance)
        {
            return FindByFeature(feature.Mass, feature.MinScanNum, feature.MaxScanNum, tolerance);
        }

        public ProteinSpectrumMatch FindByFeature(double mass, double minTime, double maxTime, Tolerance tolerance)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var minScanNum = -1;
            var maxScanNum = -1;

            for (var i = 1; i < ms1ScanNums.Length; i++)
            {
                var time = Run.GetElutionTime(ms1ScanNums[i]);
                if (minScanNum < 0 && time > minTime)
                {
                    minScanNum = ms1ScanNums[i - 1];
                }
                if (maxScanNum < 0 && time > maxTime)
                {
                    maxScanNum = ms1ScanNums[i];
                    break;
                }
            }
            return FindByFeature(mass, minScanNum, maxScanNum, tolerance);
        }

        public ProteinSpectrumMatch FindByFeature(double mass, int minScanNum, int maxScanNum, Tolerance tolerance)
        {
            var massComparerWithBinning = new MzComparerWithBinning(29); // 4ppm
            var massBin = massComparerWithBinning.GetBinNumber(mass);
            ProteinSpectrumMatch ret = null;
            for (var scanNum = minScanNum + 1; scanNum < maxScanNum; scanNum++)
            {
                if (_scanNumToMatchMap[scanNum] == null) continue;
                if (Math.Abs(massBin - massComparerWithBinning.GetBinNumber(_scanNumToMatchMap[scanNum].Mass)) <= 1)
                {
                    if (ret == null) ret = _scanNumToMatchMap[scanNum];
                    else
                    {
                        if (ret.Score < _scanNumToMatchMap[scanNum].Score) ret = _scanNumToMatchMap[scanNum];
                    }
                }
            }
            return ret;
        }

        public int CountIdentifiedScans()
        {
            return _scanNumToMatchMap.Count;
        }

        public int CountIdentifiedUniqueProteoforms()
        {
            if (ProteinSpectrumMatches == null) return 0;
            var uniqueSeq = new HashSet<string>();
            foreach (var prsm in ProteinSpectrumMatches) uniqueSeq.Add(prsm.SequenceText);
            return uniqueSeq.Count;
        }

        public int CountIdentifiedProteins()
        {
            if (ProteinSpectrumMatches == null) return 0;
            var uniqueSeq = new HashSet<string>();
            foreach (var prsm in ProteinSpectrumMatches) uniqueSeq.Add(prsm.ProteinName);
            return uniqueSeq.Count;
        }

        public readonly List<ProteinSpectrumMatch> ProteinSpectrumMatches;
        public readonly LcMsRun Run;
        public int DataId { get; private set; }
        public string DataDesc { get; private set; }
        public const double FdrCutoff = 0.01;
        private readonly Dictionary<int, ProteinSpectrumMatch> _scanNumToMatchMap;
    }
}
