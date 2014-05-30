using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public abstract class OffsetFrequencyTable
    {
        public int Charge { get; private set; }
        public int Total { get; protected set; }
        protected OffsetFrequencyTable(double searchWidth, int charge=1, double binWidth=1.005)
        {
            _offsetCounts = new Histogram<double>();
            _searchWidth = searchWidth;
            _binWidth = binWidth;
            Charge = charge;
            Total = 0;
            GenerateEdges();
        }

        public abstract void AddMatch(SpectrumMatch match);
        public abstract void AddMatches(List<SpectrumMatch> matches);

        public Probability<double>[] GetProbabilities()
        {
            var bins = _offsetCounts.Bins;
            return _offsetCounts.BinEdges.Select((t, i) => new Probability<double>(t, bins[i].Count, Total)).ToArray();
        }

        public double[] GetBinEdges()
        {
            return _offsetCounts.BinEdges;
        }

        protected double GetMinMz(double mz)
        {
            return (mz - _searchWidth);
        }

        protected double GetMaxMz(double mz)
        {
            return (mz + _searchWidth);
        }

        protected void AddOffsets(List<double> offsets)
        {
            _offsetCounts.AddData(offsets);
        }

        private void GenerateEdges()
        {
            var binEdges = new List<double>();
            for (double width = 0; width >= -1 * _searchWidth; width -= _binWidth)
            {
                binEdges.Add(width);
            }
            for (double width = 0; width < _searchWidth; width += _binWidth)
            {
                binEdges.Add(width);
            }

            binEdges = binEdges.Distinct().ToList();
            binEdges.Sort();

            _offsetCounts.BinEdges = binEdges.ToArray();
        }

        private readonly double _searchWidth;
        private readonly double _binWidth;
        private readonly Histogram<double> _offsetCounts;
    }
}
