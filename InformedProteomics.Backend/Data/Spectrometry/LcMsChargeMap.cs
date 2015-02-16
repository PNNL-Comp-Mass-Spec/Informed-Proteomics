using System;
using System.Collections;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class LcMsChargeMap
    {
        public LcMsChargeMap(LcMsRun run, Tolerance tolerance)
        {
            _run = run;
            _scanToIsolationWindow = new Dictionary<int, IsolationWindow>();

            foreach(var ms2ScanNum in _run.GetScanNumbers(2))
            {
                var isoWindow = _run.GetIsolationWindow(ms2ScanNum);
                if(isoWindow != null) _scanToIsolationWindow.Add(ms2ScanNum, isoWindow);
            }

            _tolerance = tolerance;
            _map = new Dictionary<int, BitArray>();
            _comparer = new MzComparerWithBinning(30);  // 2 ppm binning

            _sequenceMassBinToScanNumsMap = new Dictionary<int, IEnumerable<int>>();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var massBinNum = GetBinNumber(sequenceMass);
            IEnumerable<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(massBinNum, out ms2ScanNums)) return ms2ScanNums;

            return new int[0];
        }

        public void CreateMassToScanNumMap()
        {
            var minMassBin = int.MaxValue;
            var maxMassBin = int.MinValue;
            var numScans = 0;
            foreach (var entry in _map)
            {
                var massBin = entry.Key;
                var bitArray = entry.Value;

                var ms2ScanList = new List<int>();
                for (var i = 0; i < bitArray.Count; i++)
                {
                    if(bitArray[i]) ms2ScanList.Add(i + _run.MinLcScan);
                }
                if (massBin > maxMassBin) maxMassBin = massBin;
                if (massBin < minMassBin) minMassBin = massBin;
                numScans += ms2ScanList.Count;
                _sequenceMassBinToScanNumsMap.Add(massBin, ms2ScanList.ToArray());
            }
            Console.WriteLine("#MS/MS matches per sequence: {0}", numScans / (float)(maxMassBin - minMassBin + 1));

            _scanToIsolationWindow = null;
        }

        public void SetMatches(double monoIsotopicMass, int minScanNum, int maxScanNum, int minCharge, int maxCharge)
        {
            // determine bit array
            var bitArray = new BitArray(_run.MaxLcScan - _run.MinLcScan + 1);
            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                IsolationWindow isolationWindow;
                if (_scanToIsolationWindow.TryGetValue(scanNum, out isolationWindow))
                {
                    var isolationWindowTargetMz = isolationWindow.IsolationWindowTargetMz;
                    var charge = (int)Math.Round(monoIsotopicMass / isolationWindowTargetMz);
                    if (charge < minCharge || charge > maxCharge) continue;
                    var mz = Ion.GetIsotopeMz(monoIsotopicMass, charge,
                        Averagine.GetIsotopomerEnvelope(monoIsotopicMass).MostAbundantIsotopeIndex);
                    if (isolationWindow.Contains(mz))
                    {
                        bitArray.Set(scanNum-_run.MinLcScan, true);
                    }
                }
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
        private Dictionary<int, IsolationWindow> _scanToIsolationWindow;
        private readonly Dictionary<int, IEnumerable<int>> _sequenceMassBinToScanNumsMap;

        private readonly Dictionary<int, BitArray> _map;
        private readonly MzComparerWithBinning _comparer;

    }
}
