using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumScore
    {
        private readonly RankedPeaks _rankedPeaks;
        private readonly RankScore _rankScorer;
        private readonly MassErrorScore _massErrorScorer;
        private readonly Spectrum _spectrum;

        public SpectrumScore(Spectrum spectrum, string rankSetFileName, string massErrorSetFileName)
        {
            _rankedPeaks = new RankedPeaks(spectrum);
            _rankScorer = new RankScore(rankSetFileName);
            _massErrorScorer = new MassErrorScore(massErrorSetFileName);
            _spectrum = spectrum;
        }

        public double GetMassScore(IonType ionType, double mass, Tolerance tolerance)
        {
            var mz = ionType.GetMz(mass);

            // Calculate rank score
            var rank = _rankedPeaks.RankMz(mz, tolerance);
            var rankScore = _rankScorer.GetScore(ionType, rank);

            // Calculate mass error score
            var peak = _spectrum.FindPeak(mz, tolerance);
            var massErrorScore = 0.0;
            if (peak != null)
            {
                var massError = peak.Mz - mz;
                massErrorScore = _massErrorScorer.GetScore(ionType, massError);
            }
            return (rankScore + massErrorScore);
        }

        public double GetPeptideScore(string peptide, Tolerance tolerance)
        {
            var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(peptide);
            var prefixes = new List<Composition>();
            var suffixes = new List<Composition>();
            for (int i = 1; i < sequence.Count; i++)
            {
                prefixes.Add(sequence.GetComposition(0, i));
            }
            for (int i = 1; i < sequence.Count; i++)
            {
                suffixes.Add(sequence.GetComposition(i, sequence.Count));
            }
            var ionTypes = _rankScorer.IonTypes;
            var totalScore = 0.0;
            foreach (var ionType in ionTypes)
            {
                var compositions = (ionType.IsPrefixIon ? prefixes : suffixes);
                totalScore += compositions.Sum(composition => GetMassScore(ionType, composition.GetIsotopeMass(0), tolerance));
            }
            return totalScore;
        }
    }
}
