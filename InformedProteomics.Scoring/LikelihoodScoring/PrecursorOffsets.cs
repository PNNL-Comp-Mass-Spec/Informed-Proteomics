using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class PrecursorOffsets
    {
        public int Charge { get; private set; }

        private readonly Dictionary<int, List<double>> _offsets;
        private readonly double _probabilityThreshold;

        public PrecursorOffsets(int charge, double probabilityThreshold=0.0)
        {
            Charge = charge;
            _offsets = new Dictionary<int, List<double>>();
            _probabilityThreshold = probabilityThreshold;

            for (int i = 1; i <= charge; i++)
            {
                _offsets.Add(i, new List<double>());
            }
        }

        public PrecursorOffsets(List<PrecursorOffsetFrequencyTable> offsetFrequencyTables, int charge,
            double probabilityThreshold = 0.0)
        {
            Charge = charge;
            _offsets = new Dictionary<int, List<double>>();
            _probabilityThreshold = probabilityThreshold;
            for (int i = 1; i <= charge; i++)
            {
                _offsets.Add(i, new List<double>());
            }

            for (int i=0; i < offsetFrequencyTables.Count; i++)
            {
                var offsetProbabilities = offsetFrequencyTables[i].GetProbabilities();
                AddOffsets(i+1, offsetProbabilities);
            }
        }

        public void AddOffsets(int charge, IEnumerable<Probability<double>> offsetProbabilities)
        {
            if (charge > Charge)
            {
                throw new Exception("charge must be less than or equal precursor charge.");
            }

            var offsets = (from offsetProb in offsetProbabilities
                           where offsetProb.Prob >= _probabilityThreshold
                           select offsetProb.DataLabel).ToList();

            if (_offsets.ContainsKey(charge))
                _offsets[charge] = offsets;
            else
                _offsets.Add(charge, offsets);
            _offsets[charge] = _offsets[charge].Distinct().ToList(); // remove duplicates
        }

        public List<double> GetChargeOffsets(int charge)
        {
            return _offsets[charge];
        }
    }
}
