using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1BasedFilter : ISequenceFilter
    {
        public Ms1BasedFilter(
            LcMsRun run,
            int minCharge = 3, int maxCharge = 30,
            double ppmTolerance = 15)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _ppmTolerance = ppmTolerance;
            _tolerance = new Tolerance(_ppmTolerance);
            _sequenceMassBinToScanNumsMap = new Dictionary<int, IList<int>>();
            PrecomputeMostAbundantMzToMatches();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var sequenceMassBinNum = CachedLcMsRun.GetBinNumber(sequenceMass);
            IList<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(sequenceMassBinNum, out ms2ScanNums)) return ms2ScanNums;

            ms2ScanNums = new List<int>();
            var averagineEnvelope = Averagine.GetIsotopomerEnvelope(sequenceMass);
            var mostAbundantIsotopeIndex = averagineEnvelope.MostAbundantIsotopeIndex;
            for (var precursorCharge = _minCharge; precursorCharge <= _maxCharge; precursorCharge++)
            {
                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(sequenceMass, precursorCharge, mostAbundantIsotopeIndex);
                var binNumber = CachedLcMsRun.GetBinNumber(mostAbundantIsotopeMz);
                IList<ChargeAndScanNum> chargeAndScanNumList;
                if (!_mostAbundantIsotopeMzIndexToChargeAndScanNums.TryGetValue(binNumber, out chargeAndScanNumList))
                {
                    continue;
                }
                foreach (var chargeAndScanNum in chargeAndScanNumList)
                {
                    if(chargeAndScanNum.Charge == precursorCharge) ms2ScanNums.Add(chargeAndScanNum.ScanNum);
                }
            }

            _sequenceMassBinToScanNumsMap.Add(sequenceMassBinNum, ms2ScanNums);
            return ms2ScanNums;
        }

        private readonly LcMsRun _run;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly double _ppmTolerance;
        private readonly Tolerance _tolerance;

        private readonly Dictionary<int, IList<int>> _sequenceMassBinToScanNumsMap; // mass -> scan numbers
        private Dictionary<int, IList<ChargeAndScanNum>> _mostAbundantIsotopeMzIndexToChargeAndScanNums;

        private void PrecomputeMostAbundantMzToMatches()
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            _mostAbundantIsotopeMzIndexToChargeAndScanNums = new Dictionary<int, IList<ChargeAndScanNum>>();
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                var possibleSequenceMatches = GetPossibleSequenceMatches(ms2ScanNum);
                foreach (var match in possibleSequenceMatches)
                {
                    var mostAbundantIsotopeMz = match.MostAbundantIsotopeMz;
                    var charge = match.Charge;
                    var deltaMz = mostAbundantIsotopeMz * _ppmTolerance * 1e-6;
                    var minMz = mostAbundantIsotopeMz - deltaMz;
                    var maxMz = mostAbundantIsotopeMz + deltaMz;
                    var minBinNum = CachedLcMsRun.GetBinNumber(minMz);
                    var maxBinNum = CachedLcMsRun.GetBinNumber(maxMz);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        IList<ChargeAndScanNum> chargeAndScanNums;
                        if (_mostAbundantIsotopeMzIndexToChargeAndScanNums.TryGetValue(binNum,
                            out chargeAndScanNums))
                        {
                            chargeAndScanNums.Add(new ChargeAndScanNum(charge, ms2ScanNum));
                        }
                        else
                        {
                            chargeAndScanNums = new List<ChargeAndScanNum> { new ChargeAndScanNum(charge, ms2ScanNum) };
                            _mostAbundantIsotopeMzIndexToChargeAndScanNums.Add(binNum, chargeAndScanNums);
                        }
                    }
                }
            }
        }

        private IEnumerable<Ms2Match> GetPossibleSequenceMatches(int ms2ScanNumber)
        {
            var productSpec = _run.GetSpectrum(ms2ScanNumber) as ProductSpectrum;
            if (productSpec == null) return null;

            var isolationWindow = productSpec.IsolationWindow;
            var minMz = isolationWindow.MinMz;
            var maxMz = isolationWindow.MaxMz;

            var matchList = new List<Ms2Match>();
            var precursorSpec = _run.GetSpectrum(_run.GetPrecursorScanNum(ms2ScanNumber));
            if (precursorSpec != null)
            {
                var peakList = precursorSpec.GetPeakListWithin(minMz, maxMz);
                var peakListSortedByIntensity = new List<Peak>(peakList);
                peakListSortedByIntensity.Sort(new IntensityComparer());
                var remainingPeakList = new LinkedList<Peak>(peakListSortedByIntensity);
                while (remainingPeakList.Any())
                {
                    var peakMz = remainingPeakList.First.Value.Mz;
                    FindMatchedIons(peakMz, peakList, ref matchList);
                    remainingPeakList.RemoveFirst();
                }
            }

            var nextScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
            var nextSpec = _run.GetSpectrum(nextScanNum);
            if (nextSpec != null)
            {
                var peakList = nextSpec.GetPeakListWithin(minMz, maxMz);
                var peakListSortedByIntensity = new List<Peak>(peakList);
                peakListSortedByIntensity.Sort(new IntensityComparer());
                var remainingPeakList = new LinkedList<Peak>(peakListSortedByIntensity);
                while (remainingPeakList.Any())
                {
                    var peakMz = remainingPeakList.First.Value.Mz;
                    FindMatchedIons(peakMz, peakList, ref matchList);
                    remainingPeakList.RemoveFirst();
                }
            }
            return matchList;
        }

        private void FindMatchedIons(double peakMz, List<Peak> peakList, ref List<Ms2Match> matchList)
        {
            for (var charge = _maxCharge; charge >= _minCharge; charge--)
            {
                // check whether next isotope peak exists
                var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                var nextIsotopePeak = PeakListUtils.FindPeak(peakList, nextIsotopeMz, _tolerance);
                if (nextIsotopePeak == null) continue;
                matchList.Add(new Ms2Match(peakMz, charge));
            }
        }

        internal class Ms2Match
        {
            public Ms2Match(double mostAbundantIsotopeMz, int charge)
            {
                MostAbundantIsotopeMz = mostAbundantIsotopeMz;
                Charge = charge;
            }

            internal double MostAbundantIsotopeMz { get; private set; }
            internal int Charge { get; private set; }
        }

        internal class ChargeAndScanNum
        {
            public ChargeAndScanNum(int charge, int scanNum)
            {
                Charge = charge;
                ScanNum = scanNum;
            }

            internal int ScanNum { get; private set; }
            internal int Charge { get; private set; }
        }

    }
}
