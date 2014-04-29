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
            int minCharge = 3, int maxCharge = 30,
            double ppmTolerance = 10,
            double minMass = 3000.0,
            double maxMass = 50000.0,
            double isotopeCorrThreshold = 0.5,
            double chargeCorrThreshold = 0.5,
            int maxNumPeaksToConsider = 40)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _isotopeCorrThresholdThreshold = isotopeCorrThreshold;
            _chargeCorrThresholdThreshold = chargeCorrThreshold;
            MaxNumPeaksToConsider = maxNumPeaksToConsider;
            _tolerance = new Tolerance(ppmTolerance);
            _comparer = new MzComparerWithPpmTolerance(ppmTolerance);
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
        private readonly MzComparerWithPpmTolerance _comparer;
        private readonly double _chargeCorrThresholdThreshold;
        private readonly double _isotopeCorrThresholdThreshold;
        private readonly LcMsMatchMap _lcMsMatchMap;

        public void PrecomputePossibleSequenceMasses()
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                //if (ms2ScanNum != 6458) continue;
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

            //if (isolationWindow.MonoisotopicMz != null)
            //{
            //    var monoisotopicMass = isolationWindow.MonoisotopicMass;
            //    if (monoisotopicMass != null)
            //    {
            //        var xic = _run.GetExtractedIonChromatogram((double)isolationWindow.MonoisotopicMz, _tolerance, ms2ScanNumber);
            //        if (xic.Count > 2)
            //        {
            //            _lcMsMatchMap.SetMatches((double)monoisotopicMass, xic[0].ScanNum, xic[xic.Count-1].ScanNum);
            //        }
            //    }
            //}
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
            var xicThisPeak = _run.GetExtractedIonChromatogram(peakMz, _tolerance, scanNum);
            if (xicThisPeak.Count < 2) return;

            for (var charge = _maxCharge; charge >= _minCharge; charge--)
            {
                // check whether next isotope peak exists
                var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                var xicNextIsotope = _run.GetExtractedIonChromatogram(nextIsotopeMz, _tolerance, scanNum);
                if (!xicNextIsotope.Any()) continue;
                if (xicThisPeak.GetCorrelation(xicNextIsotope) < 0.7) continue;

                //var nextIsotopePeak = PeakListUtils.FindPeak(specWindow, nextIsotopeMz, _tolerance);
                //if (nextIsotopePeak == null) continue;
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

                // Charge correlation
                //var numCorrelatedCharges = 1;
                //for (var c = charge + 1; c <= charge + 10; c++)
                //{
                //    var mzChargeC = Ion.GetIsotopeMz(monoIsotopicMass, charge + 1, approxMostAbundantIsotopeIndex);
                //    var xicChargeC = _run.GetExtractedIonChromatogram(mzChargeC, _tolerance, scanNum);
                //    var corrChargeC = xicChargeC.Count >= 3 ? xicThisPeak.GetCorrelation(xicChargeC) : 0;
                //    if (corrChargeC < _chargeCorrThresholdThreshold) break;
                //    numCorrelatedCharges++;
                //}

                //for (var c = charge - 1; c >= 0; c--)
                //{
                //    var mzChargeC = Ion.GetIsotopeMz(monoIsotopicMass, charge + 1, approxMostAbundantIsotopeIndex);
                //    var xicChargeC = _run.GetExtractedIonChromatogram(mzChargeC, _tolerance, scanNum);
                //    var corrChargeC = xicChargeC.Count >= 3 ? xicThisPeak.GetCorrelation(xicChargeC) : 0;
                //    if (corrChargeC < _chargeCorrThresholdThreshold) break;
                //    numCorrelatedCharges++;
                //}

                //if (numCorrelatedCharges < 3) continue;

                var mzChargePlusOne = Ion.GetIsotopeMz(monoIsotopicMass, charge + 1, approxMostAbundantIsotopeIndex);
                var xicPlusOneCharge = _run.GetExtractedIonChromatogram(mzChargePlusOne, _tolerance, scanNum);
                var corrPlusOneCharge = xicPlusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicPlusOneCharge) : 0;

                var mzChargeMinusOne = Ion.GetIsotopeMz(monoIsotopicMass, charge - 1, approxMostAbundantIsotopeIndex);
                var xicMinusOneCharge = _run.GetExtractedIonChromatogram(mzChargeMinusOne, _tolerance, scanNum);
                var corrMinusOneCharge = xicMinusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicMinusOneCharge) : 0;

                var chargeCorr = Math.Max(corrPlusOneCharge, corrMinusOneCharge);
                if (chargeCorr < _chargeCorrThresholdThreshold) continue;                    

                _lcMsMatchMap.SetMatches(monoIsotopicMass, xicThisPeak[0].ScanNum, xicThisPeak[xicThisPeak.Count-1].ScanNum);
            }
        }
    }
}
