using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    public class CompositeScorerFactory : IFragmentScorerFactory
    {
        public CompositeScorerFactory(
            ILcMsRun run,
            IMassBinning comparer,
            AminoAcidSet aaSet,
            int minProductCharge = 1, int maxProductCharge = 20,
            double productTolerancePpm = 10,
            int isotopeOffsetTolerance = 2,
            double filteringWindowSize = 1.1
            )
            : this(run, comparer, aaSet, minProductCharge, maxProductCharge, new Tolerance(productTolerancePpm), isotopeOffsetTolerance, filteringWindowSize)
        {
        }

        public CompositeScorerFactory(
            ILcMsRun run,
            IMassBinning comparer,
            AminoAcidSet aaSet,
            int minProductCharge, int maxProductCharge,
            Tolerance productTolerance,
            int isotopeOffsetTolerance = 2,
            double filteringWindowSize = 1.1)
        {
            _run = run;
            _minProductCharge = minProductCharge;
            _maxProductCharge = maxProductCharge;
            _productTolerance = productTolerance;
            FilteringWindowSize = filteringWindowSize;
            IsotopeOffsetTolerance = isotopeOffsetTolerance;
            _ms2Scorer = new ConcurrentDictionary<int, IScorer>();
            _comparer = comparer;
            _scoringGraphFactory = new ProteinScoringGraphFactory(comparer, aaSet);
        }

        public CompositeScorerFactory(
            DPbfLcMsRun run,
            IMassBinning comparer,
            AminoAcidSet aaSet,
            int minProductCharge = 1, int maxProductCharge = 20,
            double productTolerancePpm = 10,
            int isotopeOffsetTolerance = 2,
            double filteringWindowSize = 1.1
            )
            : this(run, comparer, aaSet, minProductCharge, maxProductCharge, new Tolerance(productTolerancePpm), isotopeOffsetTolerance, filteringWindowSize)
        {
        }

        public CompositeScorerFactory(
            DPbfLcMsRun run,
            IMassBinning comparer,
            AminoAcidSet aaSet,
            int minProductCharge, int maxProductCharge,
            Tolerance productTolerance,
            int isotopeOffsetTolerance = 2,
            double filteringWindowSize = 1.1,
            PbfLcMsRun fullRun = null)
        {
            _run = run;
            _minProductCharge = minProductCharge;
            _maxProductCharge = maxProductCharge;
            _productTolerance = productTolerance;
            FilteringWindowSize = filteringWindowSize;
            IsotopeOffsetTolerance = isotopeOffsetTolerance;
            _ms2Scorer = new ConcurrentDictionary<int, IScorer>();
            _comparer = comparer;
            _scoringGraphFactory = new ProteinScoringGraphFactory(comparer, aaSet);
            _fullRun = fullRun ?? run;

            foreach (var specNum in _fullRun.AllScanNumbers.Where(x => _fullRun.GetMsLevel(x) > 1))
            {
                var spec = _fullRun.GetSpectrum(specNum) as ProductSpectrum;
                if (spec?.Peaks.Length > 0)
                {
                    var refPeakInt = CompositeScorer.GetRefIntensity(spec.Peaks);
                    _referencePeakIntensities.Add(specNum, refPeakInt);
                }
            }
        }

        private readonly IMassBinning _comparer;
        private readonly ProteinScoringGraphFactory _scoringGraphFactory;

        public double FilteringWindowSize { get; }    // 1.1
        public int IsotopeOffsetTolerance { get; }   // 2

        public IScoringGraph GetMs2ScoringGraph(int scanNum, double precursorMass)
        {
            if (precursorMass > _comparer.MaxMass || precursorMass < _comparer.MinMass)
            {
                return null;
            }

            CompositeScorerBasedOnDeconvolutedSpectrum deconvScorer;
            if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
            {
                deconvScorer = GetScorer(scanNum) as CompositeScorerBasedOnDeconvolutedSpectrum;
            }
            else
            {
                deconvScorer = GetMs2Scorer(scanNum) as CompositeScorerBasedOnDeconvolutedSpectrum;
            }

            //var deconvSpec = deconvScorer.Ms2Spectrum as DeconvolutedSpectrum;
            return _scoringGraphFactory.CreateScoringGraph(deconvScorer, precursorMass);
        }

        public IScorer GetMs2Scorer(int scanNum, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
            {
                //var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
                //var spec = _fullRun.GetSpectrum(scanNum) as ProductSpectrum;
                //if (spec == null || spec.Peaks.Length == 0)
                //    return null;
                if (!_referencePeakIntensities.TryGetValue(scanNum, out var refPeakInt))
                {
                    return null;
                }
                //var deconvolutedSpec = Deconvoluter.GetCombinedDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge, IsotopeOffsetTolerance, _productTolerance, 0.7);
                //var deconvolutedSpec = Deconvoluter.GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge,  IsotopeOffsetTolerance, FilteringWindowSize, _productTolerance);
                //return deconvolutedSpec != null ? new CompositeScorerBasedOnDeconvolutedSpectrum(deconvolutedSpec, spec, _productTolerance, _comparer, activationMethod) : null;

                return this._run.GetSpectrum(scanNum) is DeconvolutedSpectrum deconvolutedSpec ?
                    new CompositeScorerBasedOnDeconvolutedSpectrum(deconvolutedSpec, refPeakInt, _productTolerance, _comparer, activationMethod) :
                    null;
            }

            if (_ms2Scorer.TryGetValue(scanNum, out var scorer))
            {
                return scorer;
            }

            return null;
        }

        public void DeconvoluteProductSpectrum(int scanNum)
        {
            try
            {
                var scorer = GetScorer(scanNum);
                if (scorer == null)
                {
                    return;
                }

                _ms2Scorer.TryAdd(scanNum, scorer);
            }
            catch (Exception ex)
            {
                throw new Exception(string.Format("Error deconvoluting scan {0} in DeconvoluteProductSpectrum: {1}", scanNum, ex.Message), ex);
            }
        }

        public IScorer GetScorer(int scanNum, double precursorMass = 0.0, int precursorCharge = 1, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
            {
                if (_ms2Scorer.TryGetValue(scanNum, out var scorer))
                {
                    return scorer;
                }

                scorer = this.GetMs2Scorer(scanNum, activationMethod);
                this._ms2Scorer.TryAdd(scanNum, scorer);
                return scorer;
            }

            try
            {
                if (_run.GetSpectrum(scanNum) is not ProductSpectrum spec)
                {
                    return null;
                }

                var deconvolutedSpec = Deconvoluter.GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge, IsotopeOffsetTolerance, FilteringWindowSize, _productTolerance);
                return deconvolutedSpec != null ? new CompositeScorerBasedOnDeconvolutedSpectrum(deconvolutedSpec, spec, _productTolerance, _comparer) : null;
            }
            catch (Exception ex)
            {
                throw new Exception(string.Format("Error getting the scorer for scan {0} in GetScorer: {1}", scanNum, ex.Message), ex);
            }
        }

        public IScorer GetScorer(ProductSpectrum spectrum, double precursorMass, int precursorCharge, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            if (spectrum is DeconvolutedSpectrum deconSpec)
            {
                return new CompositeScorerBasedOnDeconvolutedSpectrum(deconSpec, deconSpec, _productTolerance, _comparer, deconSpec.ActivationMethod);
            }

            return new CompositeScorer(
                spectrum,
                this._productTolerance,
                activationMethod: activationMethod,
                minCharge: _minProductCharge,
                maxCharge: _maxProductCharge);
        }

        public void WriteToFile(string outputFilePath)
        {
            using var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create));

            writer.Write(_minProductCharge);
            writer.Write(_maxProductCharge);
        }

        private readonly ConcurrentDictionary<int, IScorer> _ms2Scorer;    // scan number -> scorer
        private readonly ILcMsRun _run;
        //private readonly DPbfLcMsRun _run;
        private readonly PbfLcMsRun _fullRun;
        private readonly int _minProductCharge;
        private readonly int _maxProductCharge;
        private readonly Tolerance _productTolerance;
        private readonly Dictionary<int, double> _referencePeakIntensities = new();

        /*
              * Scoring based on deconvoluted spectrum
              */
        /*
        internal class DeconvScorer : IScorer
        {
            private readonly double _prefixOffsetMass;
            private readonly double _suffixOffsetMass;
            private readonly Dictionary<int, DeconvolutedPeak> _massBinToPeakMap;
            private readonly IMassBinning _comparer;
            internal readonly DeconvolutedSpectrum DeconvolutedProductSpectrum;
            private readonly ActivationMethod _activationMethod;
            public readonly double RefIntensity;

            internal DeconvScorer(DeconvolutedSpectrum deconvolutedSpectrum, Tolerance productTolerance, IMassBinning comparer,
                double refIntensity)
            {
                DeconvolutedProductSpectrum = deconvolutedSpectrum;
                RefIntensity = refIntensity;
                _comparer = comparer;
                _activationMethod = deconvolutedSpectrum.ActivationMethod;

                if (deconvolutedSpectrum.ActivationMethod != ActivationMethod.ETD)
                {
                    _prefixOffsetMass = BaseIonType.B.OffsetComposition.Mass;
                    _suffixOffsetMass = BaseIonType.Y.OffsetComposition.Mass;
                }
                else
                {
                    _prefixOffsetMass = BaseIonType.C.OffsetComposition.Mass;
                    _suffixOffsetMass = BaseIonType.Z.OffsetComposition.Mass;
                }

                //_ionMassChkBins = new BitArray(comparer.NumberOfBins);
                _massBinToPeakMap = new Dictionary<int, DeconvolutedPeak>();

                foreach (var p in deconvolutedSpectrum.Peaks)
                {
                    var mass = p.Mz;
                    var deltaMass = productTolerance.GetToleranceAsDa(mass, 1);
                    var minMass = mass - deltaMass;
                    var maxMass = mass + deltaMass;

                    var binNum = comparer.GetBinNumber(mass);

                    if (binNum < 0)
                    {
                        binNum = comparer.GetBinNumber(minMass);
                        if (binNum < 0) binNum = comparer.GetBinNumber(maxMass);
                    }

                    // filter out
                    if (binNum < 0) continue;

                    UpdateDeconvPeak(binNum, p as DeconvolutedPeak);
                    // going up
                    for (var nextBinNum = binNum + 1; nextBinNum < comparer.NumberOfBins; nextBinNum++)
                    {
                        var nextBinMass = comparer.GetMassStart(nextBinNum);
                        if (minMass < nextBinMass && nextBinMass < maxMass) UpdateDeconvPeak(nextBinNum, p as DeconvolutedPeak); //_ionMassChkBins[nextBinNum] = true;
                        else break;
                    }

                    // going down
                    for (var prevBinNum = binNum - 1; prevBinNum < comparer.NumberOfBins; prevBinNum--)
                    {
                        var prevBinMass = comparer.GetMassEnd(prevBinNum);
                        if (minMass < prevBinMass && prevBinMass < maxMass) UpdateDeconvPeak(prevBinNum, p as DeconvolutedPeak); //_ionMassChkBins[prevBinNum] = true;
                        else break;
                    }
                }
            }

            private void UpdateDeconvPeak(int binNum, DeconvolutedPeak newPeak)
            {
                if (newPeak == null) return;
                DeconvolutedPeak existingPeak;
                if (_massBinToPeakMap.TryGetValue(binNum, out existingPeak))
                {
                    if (existingPeak.Intensity < newPeak.Intensity) _massBinToPeakMap[binNum] = newPeak;
                }
                else
                {
                    _massBinToPeakMap[binNum] = newPeak;
                }
            }

            public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
            {
                var score = 0.0;
                var prefixMass = prefixFragmentComposition.Mass + _prefixOffsetMass;
                var prefixBin = _comparer.GetBinNumber(prefixMass);
                if (prefixBin >= 0 && prefixBin < _comparer.NumberOfBins)
                {
                    DeconvolutedPeak existingPeak;
                    if (_massBinToPeakMap.TryGetValue(prefixBin, out existingPeak))
                    {
                    }
                }

                var suffixMass = suffixFragmentComposition.Mass + _suffixOffsetMass;
                var suffixBin = _comparer.GetBinNumber(suffixMass);
                if (suffixBin >= 0 && suffixBin < _comparer.NumberOfBins)
                {
                    DeconvolutedPeak existingPeak;
                    if (_massBinToPeakMap.TryGetValue(suffixBin, out existingPeak))
                    {
                        //score += _model.GetNodeScore(_activationMethod, false, suffixMass, existingPeak, RefIntensity);
                        score += 1;
                        var intScore = (existingPeak.Intensity / RefIntensity) * 10;
                        var corrScore = (existingPeak.Mass > 1300 & existingPeak.Corr > 0.7) ? (existingPeak.Corr - 0.7) : 0;
                        var distScore = (existingPeak.Mass > 1300 & existingPeak.Dist < 0.07) ? 0.3 - 3.75 * existingPeak.Dist : 0;
                        score += intScore;
                        score += corrScore;
                        score += distScore;
                    }
                }
                return score;
            }
        }
        */
    }
}
