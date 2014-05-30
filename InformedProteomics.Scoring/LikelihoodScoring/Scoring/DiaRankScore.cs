using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.Scoring
{
    public class DiaRankScore
    {
        public DiaRankScore(string fileName)
        {
            _rankScore = new RankScore(fileName);
        }

        public double GetScore(Sequence sequence, int charge, int scan, LcMsRun lcmsRun)
        {
            var spectrum = lcmsRun.GetSpectrum(scan);
            var ionTypes = _rankScore.GetIonTypes(charge);
            var filteredSpectrum = SpectrumFilter.FilterIonPeaks(sequence, spectrum, ionTypes, _tolerance);
            var match = new SpectrumMatch(sequence, filteredSpectrum, charge);
            var mass = match.PrecursorComposition.Mass;
            _rankScore.Mass = mass;
            var score = 0.0;
            var rankedPeaks = new RankedPeaks(filteredSpectrum);
            foreach (var ionType in ionTypes)
            {
                var ions = match.GetCleavageIons(ionType);
                foreach (var ion in ions)
                {
                    var rank = rankedPeaks.RankIon(ion, _tolerance);
                    score += _rankScore.GetScore(ionType, rank, charge);
                }
            }
            return score;
        }

        private readonly Tolerance _tolerance = new Tolerance(10, ToleranceUnit.Ppm);
        private readonly RankScore _rankScore;
    }
}
