using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class LcMsChargeMap
    {
        private const int MaxNumMs2ScansPerFeature = int.MaxValue;
        public LcMsChargeMap(LcMsRun run, Tolerance tolerance, int maxNumMs2ScansPerMass = MaxNumMs2ScansPerFeature)
        {
            _run = run;
            _scanToIsolationWindow = new Dictionary<int, IsolationWindow>();
            _maxNumMs2ScansPerMass = maxNumMs2ScansPerMass;

            foreach (var ms2ScanNum in _run.GetScanNumbers(2))
            {
                var isoWindow = _run.GetIsolationWindow(ms2ScanNum);
                if (isoWindow != null) _scanToIsolationWindow.Add(ms2ScanNum, isoWindow);
            }

            _tolerance = tolerance;
            _map = new Dictionary<int, BitArray>();
            _comparer = new MzComparerWithBinning(30);  // 2 ppm binning

            _sequenceMassBinToScanNumsMap = new Dictionary<int, IEnumerable<int>>();
            _scanNumToMassBin = new Dictionary<int, List<int>>();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var massBinNum = GetBinNumber(sequenceMass);
            IEnumerable<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(massBinNum, out ms2ScanNums)) return ms2ScanNums;

            return new int[0];
        }

        public IEnumerable<double> GetMatchingMass(int ms2ScanNum)
        {
            List<int> massBinNums;
            if (_scanNumToMassBin.TryGetValue(ms2ScanNum, out massBinNums)) return massBinNums.Select(s => _comparer.GetMzAverage(s));
            return new double[0];
        }

        public void CreateMassToScanNumMap()
        {
            /*
            var minMassBin = int.MaxValue;
            var maxMassBin = int.MinValue;
            var numScans = 0;
            */
            foreach (var entry in _map)
            {
                var massBin = entry.Key;
                var bitArray = entry.Value;
                var ms2ScanList = new List<int>();
                for (var i = 0; i < bitArray.Count; i++)
                {
                    if (bitArray[i])
                    {
                        var ms2ScanNum = i + _run.MinLcScan;
                        ms2ScanList.Add(ms2ScanNum);
                        List<int> massBinNums;
                        if (_scanNumToMassBin.TryGetValue(ms2ScanNum, out massBinNums))
                        {
                            massBinNums.Add(massBin);
                        }
                        else
                        {
                            _scanNumToMassBin.Add(ms2ScanNum, new List<int>() { massBin });
                        }
                    }
                }
                _sequenceMassBinToScanNumsMap.Add(massBin, ms2ScanList.ToArray());
            }
            // Console.WriteLine("#MS/MS matches per sequence: {0}", numScans / (float)(maxMassBin - minMassBin + 1));
            _scanToIsolationWindow = null;
        }

        public void SetMatches(double monoIsotopicMass, int minScanNum, int maxScanNum, int repScanNum, int minCharge, int maxCharge)
        {
            if (minScanNum < _run.MinLcScan) minScanNum = _run.MinLcScan;
            if (maxScanNum > _run.MaxLcScan) maxScanNum = _run.MaxLcScan;
            if (repScanNum < minScanNum && repScanNum > maxScanNum) return;

            // Keys are elution time, values are scan number
            var registeredMs2Scans = new List<KeyValuePair<double, int>>();

            var repRt = _run.GetElutionTime(repScanNum);
            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                IsolationWindow isolationWindow;
                if (_scanToIsolationWindow.TryGetValue(scanNum, out isolationWindow))
                {
                    var isolationWindowTargetMz = isolationWindow.IsolationWindowTargetMz;
                    var charge = (int)Math.Round(monoIsotopicMass / isolationWindowTargetMz);

                    //if (charge < minCharge || charge > maxCharge) continue;

                    var mz = Ion.GetIsotopeMz(monoIsotopicMass, charge,
                        Averagine.GetIsotopomerEnvelope(monoIsotopicMass).MostAbundantIsotopeIndex);
                    if (isolationWindow.Contains(mz))
                    {
                        var rt = _run.GetElutionTime(scanNum);
                        registeredMs2Scans.Add(new KeyValuePair<double, int>(Math.Abs(rt - repRt), scanNum));
                    }
                }
            }

            // determine bit array
            var bitArray = new BitArray(_run.MaxLcScan - _run.MinLcScan + 1);
            foreach (var e in registeredMs2Scans.OrderBy(x => x.Key).Take(_maxNumMs2ScansPerMass))
            {
                var scanNum = e.Value;
                bitArray.Set(scanNum - _run.MinLcScan, true);
            }

            var deltaMass = _tolerance.GetToleranceAsDa(monoIsotopicMass, 1);

            var minBinNum = GetBinNumber(monoIsotopicMass - deltaMass);
            var maxBinNum = GetBinNumber(monoIsotopicMass + deltaMass);

            for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
            {
                BitArray scanBitArray;
                if (!_map.TryGetValue(binNum, out scanBitArray))
                {
                    _map.Add(binNum, bitArray);
                }
                else
                {
                    scanBitArray.Or(bitArray);
                }
            }
        }

        private int GetBinNumber(double mass)
        {
            //return (int) Math.Round(mass*Constants.RescalingConstantHighPrecision);
            return _comparer.GetBinNumber(mass);
        }

        private readonly LcMsRun _run;
        private readonly Tolerance _tolerance;
        private readonly int _maxNumMs2ScansPerMass;

        private Dictionary<int, IsolationWindow> _scanToIsolationWindow;
        private readonly Dictionary<int, IEnumerable<int>> _sequenceMassBinToScanNumsMap;

        private readonly Dictionary<int, List<int>> _scanNumToMassBin;

        private readonly Dictionary<int, BitArray> _map;
        private readonly MzComparerWithBinning _comparer;
    }
}
