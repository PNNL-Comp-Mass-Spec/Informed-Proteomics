using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class PrecursorOffsetFrequencyTable
    {
        public int Charge { get; }
        public int Total { get; protected set; }
        public PrecursorOffsetFrequencyTable(double searchWidth, int charge = 1, double binWidth = 1.005)
        {
            _offsetCounts = new Histogram<double>();
            _searchWidth = searchWidth;
            _binWidth = binWidth;
            Charge = charge;
            Total = 0;
            GenerateEdges();
            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, charge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        public void AddMatch(SpectrumMatch match)
        {
            Total++;
            var ion = _precursorIonTypes[Charge - 1].GetIon(match.PeptideComposition);
            var monoIsotopicMz = ion.GetMonoIsotopicMz();
            var min = GetMinMz(monoIsotopicMz);
            var max = GetMaxMz(monoIsotopicMz);

            var peaks = match.Spectrum.Peaks;
            var offsetMzCollection = (from peak in peaks
                                      where peak.Mz >= min && peak.Mz <= max
                                      select peak.Mz - monoIsotopicMz).ToList();
            AddOffsets(offsetMzCollection);
        }

        public void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                AddMatch(match);
            }
        }

        public Probability<double>[] GetProbabilities()
        {
            var bins = _offsetCounts.Bins;
            return _offsetCounts.BinEdges.Select((t, i) => new Probability<double>(t, bins[i].Count, Total)).ToArray();
        }

        private double GetMinMz(double mz)
        {
            return (mz - _searchWidth);
        }

        private double GetMaxMz(double mz)
        {
            return (mz + _searchWidth);
        }

        private void AddOffsets(List<double> offsets)
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
        private readonly List<IonType> _precursorIonTypes;
        private readonly double _searchWidth;
        private readonly double _binWidth;
        private readonly Histogram<double> _offsetCounts;
    }
}
