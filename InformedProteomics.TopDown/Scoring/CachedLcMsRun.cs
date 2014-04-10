using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class CachedLcMsRun
    {
        public CachedLcMsRun(
            LcMsRun run, 
            int minPrecursorCharge, int maxPrecursorCharge,
            int minProductCharge, int maxProductCharge,
            double minPrecursorMz, double maxPrecursorMz, 
            Tolerance precursorTolerance,
            Tolerance productTolerance)
        {
            _run = run;
            _minPrecursorCharge = minPrecursorCharge;
            _maxPrecursorCharge = maxPrecursorCharge;
            _minProductCharge = minProductCharge;
            _maxProductCharge = maxProductCharge;
            _precursorTolerance = precursorTolerance;
            _productTolerance = productTolerance;
            _minBinNum = GetBinNumber(minPrecursorMz);
            _maxBinNum = GetBinNumber(maxPrecursorMz);

            // binNum (most abundant isotop mz), precursor charge => ms2 scan numbers
            _cachedMatchedMs2ScanNums = new List<int>[_maxBinNum - _minBinNum + 1, maxPrecursorCharge - minPrecursorCharge + 1];
            _deconvolutedMs2Spectra = new Dictionary<int, ProductSpectrum>();
            CachePrecursorMatches();
            CacheProductSpectra();
        }

        public List<Match> GetMs2Matches(Composition protCompositionWithH2O)
        {
            var matches = new List<Match>();
            for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
            {
                var mostAbundantIsotopeMz = new Ion(protCompositionWithH2O, precursorCharge).GetMostAbundantIsotopeMz();
                var binNumber = GetBinNumber(mostAbundantIsotopeMz);
                if (binNumber < _minBinNum || binNumber > _maxBinNum) continue;
                foreach (
                    var ms2ScanNum in
                        _cachedMatchedMs2ScanNums[binNumber - _minBinNum, precursorCharge - _minPrecursorCharge])
                {
                    matches.Add(new Match(precursorCharge, ms2ScanNum));
                }
            }
            return matches;
        }

        public IScorer GetMs2Scorer(int scanNum)
        {
            ProductSpectrum productSpectrum;
            if(_deconvolutedMs2Spectra.TryGetValue(scanNum, out productSpectrum)) return new DeconvScorer(productSpectrum, _productTolerance);
            return null;
        }

        private void CachePrecursorMatches()
        {
            var numBins = (float)(_maxBinNum - _minBinNum + 1);

            var cachedXic = new Dictionary<int, Xic>();
            for (var binNum = _minBinNum; binNum <= _maxBinNum; binNum++)
            {
                var mz = GetMz(binNum);
                cachedXic[binNum] = _run.GetExtractedIonChromatogram(mz, _precursorTolerance);
            }

            for (var binNum = _minBinNum; binNum <= _maxBinNum; binNum++)
            {
                if (binNum % 10000 == 0) Console.WriteLine("Caching Ms1 features {0}% done.", (binNum - _minBinNum) / numBins * 100);
                var mostAbundantIsotopeMz = GetMz(binNum);
                var ms2ScanNumArr = _run.GetFragmentationSpectraScanNums(mostAbundantIsotopeMz);
                var xicMostAbundantIsotope = cachedXic[binNum];

                var binIndex = binNum - _minBinNum;
                for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
                {
                    var matchedMs2ScanNums = new List<int>();
                    var deltaMz = Constants.C13MinusC12/precursorCharge;
                    var isotopeMinusOneBinNum = GetBinNumber(mostAbundantIsotopeMz - deltaMz);
                    var xicMostAbundantIsotopeMinusOne = isotopeMinusOneBinNum >= _minBinNum ? cachedXic[isotopeMinusOneBinNum] : null;
                    var isotopePlusOneBinNum = GetBinNumber(mostAbundantIsotopeMz + deltaMz);
                    var xicMostAbundantIsotopePlusOne = isotopePlusOneBinNum <= _maxBinNum ? cachedXic[isotopePlusOneBinNum] : null;

                    foreach (var ms2ScanNumber in ms2ScanNumArr)
                    {
                        var precursorMs1ScanNum = _run.GetPrecursorScanNum(ms2ScanNumber);
                        if (!xicMostAbundantIsotope.ContainsScanNum(precursorMs1ScanNum)) continue;
                        if (xicMostAbundantIsotopeMinusOne != null &&
                            !xicMostAbundantIsotopeMinusOne.ContainsScanNum(precursorMs1ScanNum)) continue;
                        if (xicMostAbundantIsotopePlusOne != null &&
                            !xicMostAbundantIsotopePlusOne.ContainsScanNum(precursorMs1ScanNum)) continue;
                        matchedMs2ScanNums.Add( ms2ScanNumber);
                    }
                    _cachedMatchedMs2ScanNums[binIndex, precursorCharge - _minPrecursorCharge] = matchedMs2ScanNums;
                }
            }
        }

        public void CacheProductSpectra()
        {
            //var numMs2Spectra = 0;
            //var numMonoIsotopePeaks = 0;

            for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
            {
                if (_run.GetMsLevel(scanNum) != 2) continue;

                //++numMs2Spectra;
                var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
                if (spec == null) continue;
                var peaks = spec.Peaks;

                var monoIsotopePeakList = new List<Peak>();
                var massHash = new HashSet<int>();
                for (var peakIndex = 0; peakIndex < peaks.Length; peakIndex++)
                {
                    var peak = peaks[peakIndex];

                    // Check whether peak has the maximum intensity within the window
                    var isBest = true;

                    var prevIndex = peakIndex - 1;
                    while (prevIndex >= 0)
                    {
                        var prevPeak = peaks[prevIndex];
                        if ((peak.Mz - prevPeak.Mz) > FilteringWindowSize) break;
                        if (prevPeak.Intensity > peak.Intensity)
                        {
                            isBest = false;
                            break;
                        }
                        prevIndex--;
                    }

                    if (!isBest) continue;

                    var nextIndex = peakIndex + 1;
                    while (nextIndex < peaks.Length)
                    {
                        var nextPeak = peaks[nextIndex];
                        if ((nextPeak.Mz - peak.Mz) > FilteringWindowSize) break;
                        if (nextPeak.Intensity > peak.Intensity)
                        {
                            isBest = false;
                            break;
                        }
                        nextIndex++;
                    }

                    if (!isBest) continue;

                    // peak has the maximum intensity, window = [prevIndex+1,nextIndex-1]

                    var window = new Peak[nextIndex - prevIndex - 1];
                    Array.Copy(peaks, prevIndex+1, window, 0, window.Length);
                    var windowSpectrum = new Spectrum(window, spec.ScanNum);
                    var peakMz = peak.Mz;
                    for (var productCharge = _minProductCharge; productCharge <= _maxProductCharge; productCharge++)
                    {
                        var mass = peak.Mz * productCharge;
                        var mostAbundantIsotopeIndex = Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex;

                        for (var isotopeIndex = mostAbundantIsotopeIndex - IsotopeOffsetTolerance; isotopeIndex <= mostAbundantIsotopeIndex+IsotopeOffsetTolerance; isotopeIndex++)
                        {
                            var monoIsotopeMass = (peakMz - Constants.Proton)*productCharge -
                                           isotopeIndex*Constants.C13MinusC12;
                            var massBinNum = (int) Math.Round(monoIsotopeMass*RescalingConstantHighPrecision);
                            if (massHash.Contains(massBinNum)) continue;

                            var isotopomerEnvelope = Averagine.GetIsotopomerEnvelope(monoIsotopeMass);
                            var observedPeaks = windowSpectrum.GetAllIsotopePeaks(monoIsotopeMass, productCharge, isotopomerEnvelope,
                                _productTolerance, 0.1);
                            if (observedPeaks == null) continue;

                            var envelop = isotopomerEnvelope.Envolope; 
                            var observedIntensities = new double[observedPeaks.Length];
                            
                            for (var i = 0; i < observedPeaks.Length; i++)
                            {
                                var observedPeak = observedPeaks[i];
                                observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
                            }
                            var corr = FitScoreCalculator.GetPearsonCorrelation(envelop, observedIntensities);
                            if (corr <= 0.7) continue;

                            // monoIsotopeMass is valid
                            monoIsotopePeakList.Add(new Peak(monoIsotopeMass, 1.0));
                            massHash.Add(massBinNum);
                        }
                    }
                }
                //numMonoIsotopePeaks += monoIsotopePeakList.Count;
                _deconvolutedMs2Spectra[scanNum] = new ProductSpectrum(monoIsotopePeakList, scanNum)
                {
                    ActivationMethod = spec.ActivationMethod
                };
            }
        }

        internal class DeconvScorer : IScorer
        {
            private readonly Spectrum _deconvolutedSpectrum;
            private readonly Tolerance _productTolerance;
            private readonly BaseIonType[] _baseIonTypes;
            internal DeconvScorer(ProductSpectrum deconvolutedSpectrum, Tolerance productTolerance)
            {
                _deconvolutedSpectrum = deconvolutedSpectrum;
                _productTolerance = productTolerance;
                _baseIonTypes = deconvolutedSpectrum.ActivationMethod != ActivationMethod.ETD
                    ? CorrMatchedPeakCounter.BaseIonTypesCID : CorrMatchedPeakCounter.BaseIonTypesETD;
            }

            public double GetPrecursorIonScore(Ion precursorIon)
            {
                return 0.0;
            }

            public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
            {
                var score = 0.0;

                foreach (var baseIonType in _baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                                  ? prefixFragmentComposition + baseIonType.OffsetComposition
                                  : suffixFragmentComposition + baseIonType.OffsetComposition;
                    var matchedPeak = _deconvolutedSpectrum.FindPeak(fragmentComposition.Mass, _productTolerance);
                    if (matchedPeak != null) score += matchedPeak.Intensity;
                }
                return score;
            }
        }

        private int GetBinNumber(double mass)
        {
            return (int) Math.Round(mass*RescalingConstantHighPrecision);
        }

        private double GetMz(int binNum)
        {
            return binNum/RescalingConstantHighPrecision;
        }

        private readonly List<int>[,] _cachedMatchedMs2ScanNums;
        private readonly Dictionary<int, ProductSpectrum> _deconvolutedMs2Spectra;  // scan number -> list of mono masses

        private readonly LcMsRun _run;
        private readonly int _minPrecursorCharge;
        private readonly int _maxPrecursorCharge;
        private readonly int _minProductCharge;
        private readonly int _maxProductCharge;
        private readonly int _minBinNum;
        private readonly int _maxBinNum;
        private readonly Tolerance _precursorTolerance;
        private readonly Tolerance _productTolerance;
        private const double RescalingConstantHighPrecision = Constants.RescalingConstantHighPrecision;
        private const double FilteringWindowSize = 1.1;
        private const int IsotopeOffsetTolerance = 2;

        public class Match
        {
            public Match(int precursorCharge, int scanNum)
            {
                PrecursorCharge = precursorCharge;
                ScanNum = scanNum;
            }
            public int PrecursorCharge { get; private set; }
            public int ScanNum { get; private set; }
        }
    }

}
