using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class ProductScorerBasedOnDeconvolutedSpectra
    {
        public ProductScorerBasedOnDeconvolutedSpectra(
            ILcMsRun run,
            int minProductCharge = 1, int maxProductCharge = 20,
            double productTolerancePpm = 10,
            int isotopeOffsetTolerance = 2,
            double filteringWindowSize = 1.1
            )
            : this(run, minProductCharge, maxProductCharge, new Tolerance(productTolerancePpm), isotopeOffsetTolerance, filteringWindowSize)
        {
        }

        public ProductScorerBasedOnDeconvolutedSpectra(
            ILcMsRun run,
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
        }

        public double FilteringWindowSize { get; }  // 1.1
        public int IsotopeOffsetTolerance { get; }  // 2

        public IScorer GetMs2Scorer(int scanNum)
        {
            if (_ms2Scorer.TryGetValue(scanNum, out var scorer))
            {
                return scorer;
            }

            scorer = GetScorer(scanNum);
            if (scorer == null)
            {
                return null;
            }

            return _ms2Scorer[scanNum] = scorer;
        }

        public void DeconvoluteAllProductSpectra()
        {
            foreach (var scanNum in _run.GetScanNumbers(2))
            {
                GetScorer(scanNum);
            }
        }

        public IScorer GetScorer(int scanNum)
        {
            if (!(_run.GetSpectrum(scanNum) is ProductSpectrum spec))
            {
                return null;
            }

            if (GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge, _productTolerance, CorrScoreThresholdMs2) is ProductSpectrum deconvolutedSpec)
            {
                return new DeconvScorer(deconvolutedSpec, _productTolerance);
            }

            return null;
        }

        public Spectrum GetDeconvolutedSpectrum(Spectrum spec, int minCharge, int maxCharge, Tolerance tolerance, double corrThreshold)
        {
            return GetDeconvolutedSpectrum(spec, minCharge, maxCharge, tolerance, corrThreshold, IsotopeOffsetTolerance, FilteringWindowSize);
        }

        public static Spectrum GetDeconvolutedSpectrum(
            Spectrum spec, int minCharge, int maxCharge, Tolerance tolerance, double corrThreshold,
            int isotopeOffsetTolerance, double filteringWindowSize = 1.1)
        {
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spec.ScanNum, spec.Peaks, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize, tolerance, corrThreshold);
            var peakList = new List<Peak>();
            var binHash = new HashSet<int>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = GetBinNumber(mass);
                if (!binHash.Add(binNum))
                {
                    continue;
                }

                peakList.Add(new Peak(mass, deconvolutedPeak.Intensity));
            }

            if (spec is ProductSpectrum productSpec)
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

        /*
         * Scoring based on deconvoluted spectrum
         */
        internal class DeconvScorer : IScorer
        {
            private readonly double _prefixOffsetMass;
            private readonly double _suffixOffsetMass;
            private readonly HashSet<int> _ionMassBins;
            internal DeconvScorer(ProductSpectrum deconvolutedSpectrum, Tolerance productTolerance)
            {
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

            public double GetPrecursorIonScore()
            {
                return 0.0;
            }

            public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
                AminoAcid nTerminalResidue = null,
                AminoAcid cTerminalResidue = null)
            {
                var score = 0.0;

                var prefixMass = prefixFragmentComposition.Mass + _prefixOffsetMass;
                if (_ionMassBins.Contains(GetBinNumber(prefixMass)))
                {
                    score++;
                }

                var suffixMass = suffixFragmentComposition.Mass + _suffixOffsetMass;
                if (_ionMassBins.Contains(GetBinNumber(suffixMass)))
                {
                    score++;
                }

                return score;
            }
        }

        public static int GetBinNumber(double mass)
        {
            return (int)Math.Round(mass * RescalingConstantHighPrecision);
            //return Comparer.GetBinNumber(mass);
        }

        public static double GetMz(int binNum)
        {
            return binNum / RescalingConstantHighPrecision;
            //return Comparer.GetMzAverage(binNum);
        }

        public void WriteToFile(string outputFilePath)
        {
            using var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create));

            writer.Write(_minProductCharge);
            writer.Write(_maxProductCharge);
        }

        private readonly Dictionary<int, IScorer> _ms2Scorer;    // scan number -> scorer

        private readonly ILcMsRun _run;
        private readonly int _minProductCharge;
        private readonly int _maxProductCharge;
        private readonly Tolerance _productTolerance;
        private const double RescalingConstantHighPrecision = Constants.RescalingConstantHighPrecision;
        private const double CorrScoreThresholdMs2 = 0.7;
        //private static readonly MzComparerWithBinning Comparer = new MzComparerWithBinning(29); // max error: 4ppm
    }
}
