using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class OffsetFrequencyTable
    {
        private readonly Histogram<double> _offsetCounts;

        private readonly List<IonType> _precursorIonTypes;

        private readonly int _searchWidth;
        private readonly double _binWidth;

        public int Charge { get; private set; }

        public int Total { get; private set; }

        public OffsetFrequencyTable(int searchWidth, int charge=1, double binWidth=1.005)
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

        private void GenerateEdges()
        {
            var binEdges = new List<double>();
            for (double width = 0; width >= -1*_searchWidth; width-=_binWidth)
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

        public List<OffsetProbability> OffsetFrequencies
        {
            get
            {
                var bins = _offsetCounts.Bins;
                return _offsetCounts.BinEdges.Select((t, i) => new OffsetProbability(t, bins[i].Count, Total)).ToList();
            }
        }

        public void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                Total++;
                var ion = _precursorIonTypes[Charge-1].GetIon(match.PeptideComposition);
                var monoIsotopicMz = ion.GetMonoIsotopicMz();
                var min = monoIsotopicMz - _searchWidth;
                var max = monoIsotopicMz + _searchWidth;

                var peaks = match.Spectrum.Peaks;
                var offsetMzCollection = new List<double>();

                foreach (var peak in peaks)
                {
                    if (peak.Mz >= min && peak.Mz <= max)
                    {
                        var offset = peak.Mz - monoIsotopicMz;
                        offsetMzCollection.Add(offset);
                    }
                }

                _offsetCounts.AddData(offsetMzCollection);
            }
        }
    }
}
