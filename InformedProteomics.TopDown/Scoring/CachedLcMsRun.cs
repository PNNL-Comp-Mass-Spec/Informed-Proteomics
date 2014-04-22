using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class CachedLcMsRun: ISequenceFilter
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
            _minBinNum = GetBinNumber(minPrecursorMz);
            _maxBinNum = GetBinNumber(maxPrecursorMz);
            _precursorTolerance = precursorTolerance;
            _productTolerance = productTolerance;
            _sequenceMassBinToScanNumsMap = new Dictionary<int, IList<int>>();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var sequenceMassBinNum = GetBinNumber(sequenceMass);
            IList<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(sequenceMassBinNum, out ms2ScanNums)) return ms2ScanNums;

            ms2ScanNums = new List<int>();
            var averagineEnvelope = Averagine.GetIsotopomerEnvelope(sequenceMass);
            var mostAbundantIsotopeIndex = averagineEnvelope.MostAbundantIsotopeIndex;
            for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
            {
                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(sequenceMass, precursorCharge,mostAbundantIsotopeIndex);
                var binNumber = GetBinNumber(mostAbundantIsotopeMz);
                if (binNumber < _minBinNum || binNumber > _maxBinNum) continue;
                foreach (
                    var ms2ScanNum in
                        _cachedMatchedMs2ScanNums[binNumber - _minBinNum, precursorCharge - _minPrecursorCharge])
                {
                    ms2ScanNums.Add(ms2ScanNum);
                }
            }

            _sequenceMassBinToScanNumsMap.Add(sequenceMassBinNum, ms2ScanNums);
            return ms2ScanNums;
        }

        public IScorer GetMs2Scorer(int scanNum)
        {
            //ProductSpectrum deconvolutedMs2Spectrum;
            //if(_deconvolutedMs2Spectra.TryGetValue(scanNum, out deconvolutedMs2Spectrum)) return new DeconvScorer(deconvolutedMs2Spectrum, _productTolerance);
            //return null;
            IScorer scorer;
            if (_ms2Scorer.TryGetValue(scanNum, out scorer)) return scorer;
            return null;
        }

        public void CachePrecursorMatchesBinCentric()
        {
            _cachedMatchedMs2ScanNums = new List<int>[_maxBinNum - _minBinNum + 1, _maxPrecursorCharge - _minPrecursorCharge + 1];

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
                    var deltaMz = Constants.C13MinusC12 / precursorCharge;
                    var isotopeMinusOneBinNum = GetBinNumber(mostAbundantIsotopeMz - deltaMz);
                    var xicMostAbundantIsotopeMinusOne = isotopeMinusOneBinNum >= _minBinNum ? cachedXic[isotopeMinusOneBinNum] : null;
                    var isotopePlusOneBinNum = GetBinNumber(mostAbundantIsotopeMz + deltaMz);
                    var xicMostAbundantIsotopePlusOne = isotopePlusOneBinNum <= _maxBinNum ? cachedXic[isotopePlusOneBinNum] : null;

                    foreach (var ms2ScanNumber in ms2ScanNumArr)
                    {
                        var isPrecursorValid = true;
                        var precursorMs1ScanNum = _run.GetPrecursorScanNum(ms2ScanNumber);
                        if (!xicMostAbundantIsotope.ContainsScanNum(precursorMs1ScanNum)) isPrecursorValid = false;
                        if (isPrecursorValid && xicMostAbundantIsotopeMinusOne != null &&
                            !xicMostAbundantIsotopeMinusOne.ContainsScanNum(precursorMs1ScanNum)) isPrecursorValid = false; 
                        if (isPrecursorValid && xicMostAbundantIsotopePlusOne != null &&
                            !xicMostAbundantIsotopePlusOne.ContainsScanNum(precursorMs1ScanNum)) isPrecursorValid = false;

                        if (isPrecursorValid)
                        {
                            matchedMs2ScanNums.Add(ms2ScanNumber);
                            continue;
                        }
                        var nextMs1ScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
                        if (!xicMostAbundantIsotope.ContainsScanNum(nextMs1ScanNum)) continue;
                        if (xicMostAbundantIsotopeMinusOne != null &&
                            !xicMostAbundantIsotopeMinusOne.ContainsScanNum(nextMs1ScanNum)) continue;
                        if (xicMostAbundantIsotopePlusOne != null &&
                            !xicMostAbundantIsotopePlusOne.ContainsScanNum(nextMs1ScanNum)) continue;
                        matchedMs2ScanNums.Add(ms2ScanNumber);
                    }
                    _cachedMatchedMs2ScanNums[binIndex, precursorCharge - _minPrecursorCharge] = matchedMs2ScanNums;
                }
            }
        }

        public void DeconvoluteProductSpectra()
        {
            //_deconvolutedMs2Spectra = new Dictionary<int, ProductSpectrum>();
            _ms2Scorer = new Dictionary<int, IScorer>();
            foreach (var scanNum in _run.GetScanNumbers(2))
            {
                var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
                if (spec == null) continue;

                var deconvolutedSpec = GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge, _productTolerance, CorrScoreThresholdMs2) as ProductSpectrum;
                if (deconvolutedSpec != null) _ms2Scorer[scanNum] = new DeconvScorer(deconvolutedSpec, _productTolerance);
            }
        }

        public static Spectrum GetDeconvolutedSpectrum(Spectrum spec, int minCharge, int maxCharge, Tolerance tolerance, double corrThreshold)
        {
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spec, minCharge, maxCharge, IsotopeOffsetTolerance, FilteringWindowSize, tolerance, corrThreshold);
            var peakList = new List<Peak>();
            var binHash = new HashSet<int>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = GetBinNumber(mass);
                if (!binHash.Add(binNum)) continue;
                peakList.Add(new Peak(mass, deconvolutedPeak.Intensity));
            }

            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)
            {
                return new ProductSpectrum(peakList, spec.ScanNum)
                {
                    MsLevel = spec.MsLevel,
                    ActivationMethod = productSpec.ActivationMethod,
                    IsolationWindow = productSpec.IsolationWindow
                };
            }

            return new Spectrum(peakList, spec.ScanNum);
        }

        internal class DeconvScorer : IScorer
        {
//            private readonly Spectrum _deconvolutedSpectrum;
            private readonly Tolerance _productTolerance;
            private readonly BaseIonType[] _baseIonTypes;
            private readonly HashSet<int> _ionMassBins;
            internal DeconvScorer(ProductSpectrum deconvolutedSpectrum, Tolerance productTolerance)
            {
//                _deconvolutedSpectrum = deconvolutedSpectrum;
                _productTolerance = productTolerance;
                _baseIonTypes = deconvolutedSpectrum.ActivationMethod != ActivationMethod.ETD
                    ? CorrMatchedPeakCounter.BaseIonTypesCID : CorrMatchedPeakCounter.BaseIonTypesETD;
                _ionMassBins = new HashSet<int>();
                foreach (var p in deconvolutedSpectrum.Peaks)
                {
                    var mass = p.Mz;
                    var deltaMass = productTolerance.GetToleranceAsDa(mass, 1);
                    var minMass = mass - deltaMass;
                    var maxMass = mass + deltaMass;

                    var minBinNum = GetBinNumber(minMass);
                    var maxBinNum = GetBinNumber(maxMass);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        _ionMassBins.Add(binNum);
                    }
                }
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
                    if (_ionMassBins.Contains(GetBinNumber(fragmentComposition.Mass))) score += 1;
                }
                return score;
            }
        }

        public static int GetBinNumber(double mass)
        {
            return (int) Math.Round(mass*RescalingConstantHighPrecision);
        }

        public static double GetMz(int binNum)
        {
            return binNum/RescalingConstantHighPrecision;
        }

        public void WriteToFile(string outputFilePath)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                writer.Write(_minPrecursorCharge);
                writer.Write(_maxPrecursorCharge);
                writer.Write(_minProductCharge);
                writer.Write(_maxProductCharge);
            }
        }

        private Dictionary<int, IScorer> _ms2Scorer;    // scan number -> scorer
//        private Dictionary<int, ProductSpectrum> _deconvolutedMs2Spectra;  // scan number -> list of mono masses
        private List<int>[,] _cachedMatchedMs2ScanNums;

        private readonly Dictionary<int, IList<int>> _sequenceMassBinToScanNumsMap; // mass bin -> scan numbers

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
        private const double CorrScoreThresholdMs2 = 0.7;
    }
}
