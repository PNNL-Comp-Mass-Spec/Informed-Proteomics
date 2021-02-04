using System;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Scoring;

namespace InformedProteomics.Scoring.BottomUp
{
    public class ScoredSpectrum : IScorer
    {
        public ScoredSpectrum(Spectrum spec, RankScore scorer, int charge, double massWithH2O, Tolerance tolerance)
        {
            _rankedSpec = new RankedSpectrum(spec);
            _scorer = scorer;
            _charge = charge;
            _sequenceMass = massWithH2O;
            _tolerance = tolerance;
        }

        public double GetPrecursorIonScore(Ion precursorIon)
        {
            return 0;
        }

        public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
        {
            var score = 0.0;

            foreach (var ionType in _scorer.GetIonTypes(_charge, _sequenceMass))
            {
                var monoMz = ionType.IsPrefixIon
                    ? ionType.GetMz(prefixFragmentComposition)
                    : ionType.GetMz(suffixFragmentComposition);

                //_rankedSpec.RankIon(ionType, )
                var peak = _rankedSpec.FindPeak(monoMz, _tolerance);
                if (peak == null)
                {
                    score += _scorer.GetScore(ionType, -1, _charge, _sequenceMass);  // missing
                }
                else
                {
                    score += _scorer.GetScore(ionType, peak.Rank, _charge, _sequenceMass);
                }
            }
            return score;
        }

        private readonly RankedSpectrum _rankedSpec;
        private readonly RankScore _scorer;
        private readonly int _charge;
        private readonly double _sequenceMass;
        private readonly Tolerance _tolerance;
    }

    internal class RankedSpectrum
    {
        internal RankedSpectrum(Spectrum spec)
        {
            Peaks = spec.Peaks.OrderByDescending(p => p.Intensity)
                .Select((p, i) => new RankedPeak(p.Mz, p.Intensity, i))
                .OrderBy(p => p.Mz)
                .ToArray();
        }

        internal RankedPeak[] Peaks { get; set; }

        internal RankedPeak FindPeak(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;

            var index = Array.BinarySearch(Peaks, new RankedPeak((minMz + maxMz) / 2, 0, 0));
            if (index < 0)
            {
                index = ~index;
            }

            RankedPeak bestPeak = null;
            var bestIntensity = 0.0;

            // go down
            var i = index - 1;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz <= minMz)
                {
                    break;
                }

                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeak = Peaks[i];
                }
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz >= maxMz)
                {
                    break;
                }

                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeak = Peaks[i];
                }
                ++i;
            }
            return bestPeak;
        }
    }

    internal class RankedPeak : Peak
    {
        internal RankedPeak(double mz, double intensity, int rank)
            : base(mz, intensity)
        {
            Rank = rank;
        }

        internal int Rank { get; }
    }
}
