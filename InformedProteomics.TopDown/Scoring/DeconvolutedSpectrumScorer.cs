using System.Collections;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    public class DeconvolutedSpectrumScorer
    {
        public DeconvolutedSpectrumScorer(
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

        public DeconvolutedSpectrumScorer(
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
            var deconvScorer = GetMs2Scorer(scanNum) as DeconvScorer;
            return _scoringGraphFactory.CreateScoringGraph(deconvScorer.DeconvolutedProductSpectrum, precursorMass);
        }

        public IScorer GetMs2Scorer(int scanNum)
        {
            IScorer scorer;
            if (_ms2Scorer.TryGetValue(scanNum, out scorer)) return scorer;
            scorer = GetScorer(scanNum);
            if (scorer == null) return null;

            lock (_ms2Scorer)
            {
                _ms2Scorer.Add(scanNum, scorer);
            }
            return scorer;
        }

        public void DeconvoluteAllProductSpectra()
        {
            foreach (var scanNum in _run.GetScanNumbers(2))
            {
                GetMs2Scorer(scanNum);
            }
        }

        public IScorer GetScorer(int scanNum)
        {
            var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
            if (spec == null) return null;
            var deconvolutedSpec = GetDeconvolutedSpectrum(spec, _minProductCharge, _maxProductCharge, _productTolerance, CorrScoreThresholdMs2, IsotopeOffsetTolerance, FilteringWindowSize) as ProductSpectrum;
            if (deconvolutedSpec != null) return new DeconvScorer(deconvolutedSpec, _productTolerance, _comparer);
            return null;
        }

        public Spectrum GetDeconvolutedSpectrum(Spectrum spec, int minCharge, int maxCharge, Tolerance tolerance, double corrThreshold,
                                                       int isotopeOffsetTolerance, double filteringWindowSize = 1.1)
        {
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spec, minCharge, maxCharge, isotopeOffsetTolerance, filteringWindowSize, tolerance, corrThreshold);
            var peakList = new List<Peak>();
            var binHash = new HashSet<int>();
            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = _comparer.GetBinNumber(mass);
                if (binNum < 0) continue; // ignore peaks that cannot be a combination of amino acid masses

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

        /*
         * Scorign based on deconvoluted spectrum
         */
        internal class DeconvScorer : IScorer
        {
            private readonly double _prefixOffsetMass;
            private readonly double _suffixOffsetMass;
            private readonly BitArray _ionMassChkBins;
            private readonly IMassBinning _comparer;
            internal readonly ProductSpectrum DeconvolutedProductSpectrum;

            internal DeconvScorer(ProductSpectrum deconvolutedSpectrum, Tolerance productTolerance, IMassBinning comparer)
            {
                DeconvolutedProductSpectrum = deconvolutedSpectrum;
                _comparer = comparer;
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
                
                _ionMassChkBins = new BitArray(comparer.NumberOfBins);

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

                    _ionMassChkBins[binNum] = true;
                    // going up
                    for (var nextBinNum = binNum + 1; nextBinNum < comparer.NumberOfBins; nextBinNum++)
                    {
                        var nextBinMass = comparer.GetMassStart(nextBinNum);
                        if (minMass < nextBinMass && nextBinMass < maxMass) _ionMassChkBins[nextBinNum] = true;
                        else break;
                    }

                    // going down
                    for (var prevBinNum = binNum - 1; prevBinNum < comparer.NumberOfBins; prevBinNum--)
                    {
                        var prevBinMass = comparer.GetMassEnd(prevBinNum);
                        if (minMass < prevBinMass && prevBinMass < maxMass) _ionMassChkBins[prevBinNum] = true;
                        else break;
                    }

                    /*
                    var minBinNum = comparer.GetBinNumber(minMass);
                    var maxBinNum = comparer.GetBinNumber(maxMass);
                                        for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        _ionMassBins.Add(binNum);
                    }*/
                }
            }

            public double GetPrecursorIonScore(Ion precursorIon)
            {
                return 0.0;
            }

            public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
            {
                var score = 0.0;

                var prefixMass = prefixFragmentComposition.Mass + _prefixOffsetMass;
                var prefixBin = _comparer.GetBinNumber(prefixMass);
                if (prefixBin >= 0 && prefixBin < _ionMassChkBins.Length && _ionMassChkBins[prefixBin]) score += 1;
                
                var suffixMass = suffixFragmentComposition.Mass + _suffixOffsetMass;
                var suffixBin = _comparer.GetBinNumber(suffixMass);
                if (suffixBin >= 0 && suffixBin < _ionMassChkBins.Length &&_ionMassChkBins[suffixBin]) score += 1;

                return score;
            }
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
        private const double CorrScoreThresholdMs2 = 0.7;

        //private static readonly MzComparerWithBinning Comparer = new MzComparerWithBinning(29); // max error: 4ppm
    }
}
