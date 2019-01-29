using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class PrecursorOffsets
    {
        public int Charge { get; private set; }

        public PrecursorOffsets(int charge, double searchWidth=100, double probabilityThreshold=0.0, double binWidth=1.005)
        {
            Charge = charge;
            _offsets = new Dictionary<int, PrecursorOffsetFrequencyTable>();
            _probabilityThreshold = probabilityThreshold;

            for (var i = 1; i <= charge; i++)
            {
                _offsets.Add(i, new PrecursorOffsetFrequencyTable(searchWidth/(2*i), i, binWidth/i));
            }
        }

        public void AddMatch(SpectrumMatch match)
        {
            var charge = match.PrecursorCharge;
            if (!_offsets.ContainsKey(charge)) throw new Exception("Invalid charge.");
            for (var i = charge; i > 0; i--)
            {
                _offsets[i].AddMatches(new List<SpectrumMatch>{match});
            }
        }

        public void AddMatches(IEnumerable<SpectrumMatch> matches)
        {
            foreach (var match in matches)  AddMatch(match);
        }

        public List<double> GetChargeOffsets(int charge)
        {
            if (charge > Charge || charge < 1)  throw new Exception("Charge must be between 1 and "+Charge);

            var bins = _offsets[charge].GetProbabilities();
            var offsets = (from offsetProb in bins
                           where offsetProb.Prob >= _probabilityThreshold
                           select offsetProb.Label).ToList();
            return offsets.Distinct().ToList();
        }

        private readonly Dictionary<int, PrecursorOffsetFrequencyTable> _offsets;
        private readonly double _probabilityThreshold;
    }
}
