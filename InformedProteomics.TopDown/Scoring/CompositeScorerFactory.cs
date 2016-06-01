using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.TopDown.Scoring
{
    public class CompositeScorerFactory
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
            _ms2Scorer = new Dictionary<int, IScorer>();
            _comparer = comparer;
            _scoringGraphFactory = new ProteinScoringGraphFactory(comparer, aaSet);
        }

        private readonly IMassBinning _comparer;
        private readonly ProteinScoringGraphFactory _scoringGraphFactory;

        public double FilteringWindowSize { get; private set; }    // 1.1
        public int IsotopeOffsetTolerance { get; private set; }   // 2

        public IScoringGraph GetMs2ScoringGraph(int scanNum, double precursorMass)
        {
            if (precursorMass > _comparer.MaxMass || precursorMass < _comparer.MinMass) return null;

            var deconvScorer = GetMs2Scorer(scanNum) as CompositeScorerBasedOnDeconvolutedSpectrum;
            //var deconvSpec = deconvScorer.Ms2Spectrum as DeconvolutedSpectrum;
            return _scoringGraphFactory.CreateScoringGraph(deconvScorer, precursorMass);
        }

        public IScorer GetMs2Scorer(int scanNum)
        {
            IScorer scorer;
            if (_ms2Scorer.TryGetValue(scanNum, out scorer)) return scorer;
            return null;
        }

        public void DeconvonluteProductSpectrum(int scanNum)
        {
            var scorer = GetScorer(scanNum);
            if (scorer == null) return;
            lock (_ms2Scorer)
            {
                _ms2Scorer.Add(scanNum, scorer);
            }
        }

        public IScorer GetScorer(int scanNum)
        {
            var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
            if (spec == null) return null;
            var deconvolutedSpec = Deconvoluter.GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge,  IsotopeOffsetTolerance, FilteringWindowSize, _productTolerance);
            return deconvolutedSpec != null ? new CompositeScorerBasedOnDeconvolutedSpectrum(deconvolutedSpec, spec, _productTolerance, _comparer) : null;
        }

        public void WriteToFile(string outputFilePath)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                writer.Write(_minProductCharge);
                writer.Write(_maxProductCharge);
            }
        }

        private readonly Dictionary<int, IScorer> _ms2Scorer;    // scan number -> scorer
        private readonly ILcMsRun _run;
        private readonly int _minProductCharge;
        private readonly int _maxProductCharge;
        private readonly Tolerance _productTolerance;

        /*
              * Scorign based on deconvoluted spectrum
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
