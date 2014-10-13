using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1IsotopeAndChargeCorrFilter : ISequenceFilter
    {
        public Ms1IsotopeAndChargeCorrFilter(
            LcMsRun run,
            Tolerance tolerance,
            int minCharge = 3, int maxCharge = 30,
            double minMass = 3000.0,
            double maxMass = 50000.0,
            double isotopeCorrThreshold = 0.7,
            double chargeCorrThreshold = 0.7,
            double mostAbundantPlusOneIsotopeCorrThreshold = 0.7,
            int maxNumPeaksToConsider = 40)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _isotopeCorrThresholdThreshold = isotopeCorrThreshold;
            _chargeCorrThresholdThreshold = chargeCorrThreshold;
            _mostAbundantPlusOneIsotopeCorrThreshold = mostAbundantPlusOneIsotopeCorrThreshold;
            MaxNumPeaksToConsider = maxNumPeaksToConsider;
            _tolerance = tolerance;
            _comparer = new MzComparerWithTolerance(tolerance);
            _lcMsMatchMap = new LcMsMatchMap();
            PrecomputePossibleSequenceMasses();
            _lcMsMatchMap.CreateSequenceMassToMs2ScansMap(_run, _tolerance, minMass, maxMass);
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            return _lcMsMatchMap.GetMatchingMs2ScanNums(sequenceMass, _tolerance, _run);
        }

        public int MaxNumPeaksToConsider { get; private set; }

        private readonly LcMsRun _run;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly Tolerance _tolerance;
        private readonly MzComparerWithTolerance _comparer;
        private readonly double _chargeCorrThresholdThreshold;
        private readonly double _isotopeCorrThresholdThreshold;
        private readonly double _mostAbundantPlusOneIsotopeCorrThreshold;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public void PrecomputePossibleSequenceMasses()
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                //if (ms2ScanNum != 52468) continue;
                SetLcMsMatches(ms2ScanNum);
            }
        }

        private void SetLcMsMatches(int ms2ScanNumber)
        {
            var productSpec = _run.GetSpectrum(ms2ScanNumber) as ProductSpectrum;
            if (productSpec == null) return;

            var isolationWindow = productSpec.IsolationWindow;
            var minMz = isolationWindow.MinMz;
            var maxMz = isolationWindow.MaxMz;

            List<Peak> precursorPeakList = null;
            List<Peak> precursorSpecWindow = null;
            var precursorSpec = _run.GetSpectrum(_run.GetPrecursorScanNum(ms2ScanNumber));
            if (precursorSpec != null)
            {
                precursorPeakList = precursorSpec.GetPeakListWithin(minMz, maxMz);
                precursorSpecWindow = precursorSpec.GetPeakListWithin(minMz-1.0, maxMz+1.0);
            }

            List<Peak> nextMs1PeakList = null;
            List<Peak> nextMs1SpecWindow = null;
            var nextScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
            var nextSpec = _run.GetSpectrum(nextScanNum);
            if (nextSpec != null)
            {
                nextMs1PeakList = nextSpec.GetPeakListWithin(minMz, maxMz);
                nextMs1SpecWindow = nextSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
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
            else return;

            // Sort by intensity
            peakList.Sort(new IntensityComparer());
            var remainingPeakList = new LinkedList<Peak>(peakList);

            SetLcMsMatches(remainingPeakList, ms2ScanNumber, precursorSpecWindow, nextMs1SpecWindow, MaxNumPeaksToConsider);
        }

        private void SetLcMsMatches(LinkedList<Peak> remainingPeakList, int scanNum, IList<Peak> precursorSpecWindow, IList<Peak> nextMs1SpecWindow, int numPeaksToConsider)
        {
            var numPeaksConsidered = 0;
            while(remainingPeakList.Any())
            {
                var peakWithHighestIntensity = remainingPeakList.First.Value;
                var peakMz = peakWithHighestIntensity.Mz;
                SetLcMsMatches(peakMz, scanNum, precursorSpecWindow, nextMs1SpecWindow);
                if (++numPeaksConsidered >= numPeaksToConsider) break;
                remainingPeakList.RemoveFirst();
            }
        }

        private void SetLcMsMatches(double peakMz, int scanNum, IList<Peak> precursorSpecWindow, IList<Peak> nextMs1SpecWindow)
        {
            var xicThisPeak = _run.GetPrecursorExtractedIonChromatogram(peakMz, _tolerance, scanNum);
            if (xicThisPeak.Count < 2) return;

            for (var charge = _maxCharge; charge >= _minCharge; charge--)
            {
                // check whether next isotope peak exists
                var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                var xicNextIsotope = _run.GetPrecursorExtractedIonChromatogram(nextIsotopeMz, _tolerance, scanNum);
                if (!xicNextIsotope.Any()) continue;
                if (xicThisPeak.GetCorrelation(xicNextIsotope) < _mostAbundantPlusOneIsotopeCorrThreshold) continue;

                var mostAbundantIsotopeMass = (peakMz - Constants.Proton) * charge;
                var averagineIsoEnv = Averagine.GetIsotopomerEnvelope(mostAbundantIsotopeMass);
                var approxMostAbundantIsotopeIndex = averagineIsoEnv.MostAbundantIsotopeIndex;
                var monoIsotopicMass = mostAbundantIsotopeMass - approxMostAbundantIsotopeIndex * Constants.C13MinusC12;

                // Isotope correlation
                var averagineIsotopeProfile = Averagine.GetTheoreticalIsotopeProfile(monoIsotopicMass, charge);
                var precursorIsotopeCorr = precursorSpecWindow != null ? PeakListUtils.GetPearsonCorrelation(precursorSpecWindow, averagineIsotopeProfile, _comparer) : 0;
                var nextMs1IsotopeCorr = nextMs1SpecWindow != null ? PeakListUtils.GetPearsonCorrelation(nextMs1SpecWindow, averagineIsotopeProfile, _comparer) : 0;
                var isotopeCorr = Math.Max(precursorIsotopeCorr, nextMs1IsotopeCorr);
                if (isotopeCorr < _isotopeCorrThresholdThreshold) continue;

                if (_chargeCorrThresholdThreshold > 0.0)
                {
                    var mzChargePlusOne = Ion.GetIsotopeMz(monoIsotopicMass, charge + 1, approxMostAbundantIsotopeIndex);
                    var xicPlusOneCharge = _run.GetPrecursorExtractedIonChromatogram(mzChargePlusOne, _tolerance, scanNum);
                    var corrPlusOneCharge = xicPlusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicPlusOneCharge) : 0;

                    double corrMinusOneCharge;
                    if (charge > 1)
                    {
                        var mzChargeMinusOne = Ion.GetIsotopeMz(monoIsotopicMass, charge - 1, approxMostAbundantIsotopeIndex);
                        var xicMinusOneCharge = _run.GetPrecursorExtractedIonChromatogram(mzChargeMinusOne, _tolerance, scanNum);
                        corrMinusOneCharge = xicMinusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicMinusOneCharge) : 0;
                    }
                    else
                    {
                        corrMinusOneCharge = 0.0;
                    }

                    var chargeCorr = Math.Max(corrPlusOneCharge, corrMinusOneCharge);
                    if (chargeCorr < _chargeCorrThresholdThreshold) continue;                               
                }

                _lcMsMatchMap.SetMatches(monoIsotopicMass, xicThisPeak[0].ScanNum, xicThisPeak[xicThisPeak.Count-1].ScanNum);
            }
        }
    }
}
