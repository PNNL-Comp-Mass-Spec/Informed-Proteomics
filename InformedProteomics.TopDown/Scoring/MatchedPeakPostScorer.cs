using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.TopDown.Scoring
{
    public class MatchedPeakPostScorer
    {
        public MatchedPeakPostScorer(Tolerance tolerance, int minCharge, int maxCharge)
        {
            _tolerance = tolerance;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _rankingInfo = new Dictionary<int, int[]>();
        }

        public double ComputeScore(ProductSpectrum ms2Spec, Sequence sequence)
        {
            
            _baseIonTypes = ms2Spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            _sequence = sequence;
            _ms2Spec = ms2Spec;

            if (!_rankingInfo.TryGetValue(ms2Spec.ScanNum, out _peakRanking))
            {
                _peakRanking = ArrayUtil.GetRankings(ms2Spec.Peaks.Select(peak => peak.Intensity));
                _rankingInfo.Add(ms2Spec.ScanNum, _peakRanking);
            }

            FindMatchedPeaks();

            var score = GetRankSumScore();

            return score;
        }

        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly Dictionary<int, int[]> _rankingInfo;
        private BaseIonType[] _baseIonTypes;

        private const double RelativeIsotopeIntensityThreshold = 0.8;
        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        private static readonly MzComparerWithBinning Comparer;
        
        private int[] _peakRanking;
        private ProductSpectrum _ms2Spec;
        private Sequence _sequence;
        private List<int> _prefixIonPeakIndex;
        private List<int> _suffixIonPeakIndex;
        private int _nTheoreticalIonPeaks;
        private int _nObservedIonPeaks;

        private double GetRankSumScore()
        {
            var rankSum = 0d;
            var nMatchedIons = _prefixIonPeakIndex.Count + _suffixIonPeakIndex.Count;
            var nObservedPeaks = _ms2Spec.Peaks.Length;

            foreach (var peakIndex in _prefixIonPeakIndex) rankSum += _peakRanking[peakIndex];
            foreach (var peakIndex in _suffixIonPeakIndex) rankSum += _peakRanking[peakIndex];

            var pvalue = FitScoreCalculator.GetRankSumPvalue(nObservedPeaks, nMatchedIons, rankSum);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        private double GetHyperGeometricScore()
        {
            var nPossiblePeaks = Comparer.GetBinNumber(_ms2Spec.Peaks.Last().Mz) - Comparer.GetBinNumber(_ms2Spec.Peaks.First().Mz) + 1;
            var nObservedPeaks = _ms2Spec.Peaks.Length;
            var pvalue = FitScoreCalculator.GetHyperGeometricPvalue(nPossiblePeaks, nObservedPeaks, _nTheoreticalIonPeaks, _nObservedIonPeaks);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }


        private void FindMatchedPeaks()
        {
            var cleavages = _sequence.GetInternalCleavages();

            _prefixIonPeakIndex = new List<int>(); // list of cleavage indices
            _suffixIonPeakIndex = new List<int>();
            _nTheoreticalIonPeaks = 0;
            _nObservedIonPeaks = 0;

            int index = 0; // cleavage index 
            foreach (var c in cleavages)
            {
                foreach (var baseIonType in _baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                                  ? c.PrefixComposition + baseIonType.OffsetComposition
                                  : c.SuffixComposition + baseIonType.OffsetComposition;

                    for (var charge = _minCharge; charge <= _maxCharge; charge++)
                    {
                        var ion = new Ion(fragmentComposition, charge);
                        int baseIsotopePeakIndex;
                        int nIsotopes;
                        int nMatchedIsotopes;

                        if (FindIon(ion, _tolerance, RelativeIsotopeIntensityThreshold, out baseIsotopePeakIndex, out nIsotopes, out nMatchedIsotopes))
                        {
                            if (baseIonType.IsPrefix) _prefixIonPeakIndex.Add(baseIsotopePeakIndex);
                            else _suffixIonPeakIndex.Add(baseIsotopePeakIndex);
                            _nObservedIonPeaks++;
                        }
                        //_nObservedIonPeaks += nMatchedIsotopes;
                        //_nTheoreticalIonPeaks += nIsotopes;
                        _nTheoreticalIonPeaks++;
                    }
                }
                index++;
            }
        }

        private bool FindIon(Ion ion, Tolerance tolerance, double relativeIntensityThreshold, out int baseIsotopePeakIndex, out int nIsotopes, out int nMatchedIsotopes)
        {
            //matchedPeakIndex = new List<int>();
            var baseIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
            var baseIsotopMz = ion.GetIsotopeMz(baseIsotopeIndex);
            baseIsotopePeakIndex = _ms2Spec.FindPeakIndex(baseIsotopMz, tolerance);

            nIsotopes = isotopomerEnvelope.Select(x => x >= relativeIntensityThreshold).Count();
            nMatchedIsotopes = 0;

            if (baseIsotopePeakIndex < 0) return false;
            //if (baseIsotopePeakIndex < 0) baseIsotopePeakIndex = ~baseIsotopePeakIndex;
            nMatchedIsotopes++;

            // go down
            var peakIndex = baseIsotopePeakIndex;
            //matchedPeakIndex.Add(peakIndex);
            for (var isotopeIndex = baseIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;

                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex - 1; i >= 0; i--)
                {
                    var peakMz = _ms2Spec.Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        //peakIndex = i;
                        //break;
                        return false;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        //matchedPeakIndex.Add(peakIndex);
                        nMatchedIsotopes++;
                        break;
                    }
                }
            }

            // go up
            peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;

                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex + 1; i < _ms2Spec.Peaks.Length; i++)
                {
                    var peakMz = _ms2Spec.Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        //peakIndex = i;
                        //break;
                        return false;
                    }
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        //matchedPeakIndex.Add(peakIndex);
                        nMatchedIsotopes++;
                        break;
                    }
                }
            }

            return true;
        }

        static MatchedPeakPostScorer()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
            Comparer = new MzComparerWithBinning(28);
        }
    }
}
