using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1IsotopeMostAbundantPlusOneFilter : ISequenceFilter
    {
        public Ms1IsotopeMostAbundantPlusOneFilter(
            InMemoryLcMsRun run,
            int minCharge = 3, int maxCharge = 30,
            double ppmTolerance = 10,
            double minMass = 3000.0,
            double maxMass = 50000.0,
            int maxNumPeaksToConsider = 40)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            MaxNumPeaksToConsider = maxNumPeaksToConsider;
            _tolerance = new Tolerance(ppmTolerance);
            _comparer = new MzComparerWithTolerance(ppmTolerance);
            _lcMsMatchMap = new LcMsMatchMap();
            PrecomputePossibleSequenceMasses();
            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _tolerance, minMass, maxMass);
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass);
        }

        public int MaxNumPeaksToConsider { get; }

        private readonly InMemoryLcMsRun _run;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly Tolerance _tolerance;
        private readonly MzComparerWithTolerance _comparer;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public void PrecomputePossibleSequenceMasses()
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                SetLcMsMatches(ms2ScanNum);
            }
        }

        private void SetLcMsMatches(int ms2ScanNumber)
        {
            if (!(_run.GetSpectrum(ms2ScanNumber) is ProductSpectrum productSpec))
            {
                return;
            }

            var isolationWindow = productSpec.IsolationWindow;
            var minMz = isolationWindow.MinMz;
            var maxMz = isolationWindow.MaxMz;

            List<Peak> precursorPeakList = null;
            var precursorSpec = _run.GetSpectrum(_run.GetPrecursorScanNum(ms2ScanNumber));
            if (precursorSpec != null)
            {
                precursorPeakList = precursorSpec.GetPeakListWithin(minMz, maxMz);
            }

            List<Peak> nextMs1PeakList = null;
            var nextScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
            var nextSpec = _run.GetSpectrum(nextScanNum);
            if (nextSpec != null)
            {
                nextMs1PeakList = nextSpec.GetPeakListWithin(minMz, maxMz);
            }

            List<Peak> peakList;
            if (precursorPeakList != null && nextMs1PeakList != null)
            {
                peakList = PeakListUtils.Sum(precursorPeakList, nextMs1PeakList, _comparer);
            }
            else if (precursorPeakList != null)
            {
                peakList = precursorPeakList;
            }
            else if (nextMs1PeakList != null)
            {
                peakList = nextMs1PeakList;
            }
            else
            {
                return;
            }

            // Sort by intensity
            peakList.Sort(new IntensityComparer());
            var remainingPeakList = new LinkedList<Peak>(peakList);

            SetLcMsMatches(remainingPeakList, ms2ScanNumber, MaxNumPeaksToConsider);
        }

        private void SetLcMsMatches(LinkedList<Peak> remainingPeakList, int scanNum, int numPeaksToConsider)
        {
            var numPeaksConsidered = 0;
            while (remainingPeakList.Count > 0)
            {
                var peakWithHighestIntensity = remainingPeakList.First.Value;
                var peakMz = peakWithHighestIntensity.Mz;
                SetLcMsMatches(peakMz, scanNum);
                if (++numPeaksConsidered >= numPeaksToConsider)
                {
                    break;
                }

                remainingPeakList.RemoveFirst();
            }
        }

        private void SetLcMsMatches(double peakMz, int scanNum)
        {
            var xicThisPeak = _run.GetPrecursorExtractedIonChromatogram(peakMz, _tolerance, scanNum);
            if (xicThisPeak.Count < 2)
            {
                return;
            }

            for (var charge = _maxCharge; charge >= _minCharge; charge--)
            {
                // check whether next isotope peak exists
                var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                var xicNextIsotope = _run.GetPrecursorExtractedIonChromatogram(nextIsotopeMz, _tolerance, scanNum);
                if (!xicNextIsotope.Any())
                {
                    continue;
                }

                var mostAbundantIsotopeMass = (peakMz - Constants.Proton) * charge;
                var averagineIsoEnv = Averagine.GetIsotopomerEnvelope(mostAbundantIsotopeMass);
                var approxMostAbundantIsotopeIndex = averagineIsoEnv.MostAbundantIsotopeIndex;
                var monoIsotopicMass = mostAbundantIsotopeMass - approxMostAbundantIsotopeIndex * Constants.C13MinusC12;

                _lcMsMatchMap.SetMatches(monoIsotopicMass, xicThisPeak[0].ScanNum, xicThisPeak[xicThisPeak.Count - 1].ScanNum);
            }
        }
    }
}
