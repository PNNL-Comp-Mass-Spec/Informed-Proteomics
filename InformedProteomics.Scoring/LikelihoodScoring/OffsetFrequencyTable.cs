using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class OffsetFrequencyTable
    {
        private readonly Histogram<double> _offsetCounts;

        private readonly IonType _precursorIonType;

        private readonly int _searchWidth;
        private readonly double _binWidth;

        public OffsetFrequencyTable(int searchWidth, double binWidth=1.005, int maxCharge=4)
        {
            _offsetCounts = new Histogram<double>();
            _searchWidth = searchWidth;
            _binWidth = binWidth;
            GenerateEdges();

            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, maxCharge);

            _precursorIonType = ionTypeFactory.GetIonType("y");
        }

        private void GenerateEdges()
        {
            var binEdges = new List<double>();
            for (double i = 0; i >= -1*_searchWidth; i-=_binWidth)
            {
                binEdges.Add(i);
            }
            for (double i = 0; i < _searchWidth; i += _binWidth)
            {
                binEdges.Add(i);
            }

            binEdges = binEdges.Distinct().ToList();
            binEdges.Sort();

            _offsetCounts.BinEdges = binEdges.ToArray();
        }

        public List<OffsetProbability> OffsetFrequencies
        {
            get
            {
                var offsetFrequencies = new List<OffsetProbability>();
                var bins = _offsetCounts.Bins;
                for (int i = 0; i < _offsetCounts.BinEdges.Length; i++)
                {
                    if (bins[i].Count > 0)
                        offsetFrequencies.Add(new OffsetProbability(_offsetCounts.BinEdges[i], bins[i].Count));
                }
                return offsetFrequencies;
            }
        }

        public void AddMatches(List<SpectrumMatch> matches, Tolerance tolerance)
        {
            foreach (var match in matches)
            {
                var ion = _precursorIonType.GetIon(match.PeptideComposition);
                var monoIsotopicMz = ion.GetMonoIsotopicMz();
                var min = monoIsotopicMz - _searchWidth;
                var max = monoIsotopicMz + _searchWidth;
                int basePeakIndex;
                if (tolerance.GetUnit() == ToleranceUnit.Da)
                {
                    var toleranceValue = tolerance.GetValue();
                    basePeakIndex = match.Spectrum.FindPeakIndex(monoIsotopicMz - toleranceValue, monoIsotopicMz + toleranceValue);
                }
                else
                {
                    basePeakIndex = match.Spectrum.FindPeakIndex(monoIsotopicMz, tolerance);
                }
                var peaks = match.Spectrum.Peaks;
                var offsetMzCollection = new List<double>();

                if (basePeakIndex <= 0 || basePeakIndex >= peaks.Length) continue;
                
                for (int i = basePeakIndex;  i >= 0 && peaks[i].Mz >= min; i--)
                {
                    var currentMz = peaks[i].Mz;
                    var offset = currentMz - monoIsotopicMz;

                    offsetMzCollection.Add(offset);
                }

                for (int i = basePeakIndex; i < peaks.Length && peaks[i].Mz <= max; i++)
                {
                    var currentMz = peaks[i].Mz;
                    var offset = currentMz - monoIsotopicMz;

                    offsetMzCollection.Add(offset);
                }

                _offsetCounts.AddData(offsetMzCollection);
            }
        }
    }
}
