using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class PrecursorFilter
    {
        private readonly Dictionary<int, PrecursorOffsets> _offsets; 
        private readonly List<IonType> _precursorIonTypes;
        private readonly Tolerance _tolerance;
        private readonly double _searchWidth;
        private readonly int _retentionCount;

        public PrecursorFilter(Dictionary<int, PrecursorOffsets> offsets, Tolerance tolerance,
                               double searchWidth, int retentionCount)
        {
            _offsets = offsets;
            _tolerance = tolerance;
            _searchWidth = searchWidth;
            _retentionCount = retentionCount;

            var maxCharge = offsets.Keys.Max();

            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, maxCharge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        private SpectrumMatch Filter(SpectrumMatch specMatch, bool filterNoise)
        {
            var charge = specMatch.PrecursorCharge;
            var peaks = specMatch.Spectrum.Peaks.ToList();
            var offsetLists = _offsets[charge];

            for (int i = 0; i < _precursorIonTypes.Count; i++)
            {
                var ion = _precursorIonTypes[i].GetIon(specMatch.PeptideComposition);
                var mz = ion.GetMonoIsotopicMz();

                if (filterNoise)
                    peaks = (SpectrumFilter.FilterNoisePeaks(specMatch.Spectrum)).Peaks.ToList();

                var offsets = offsetLists.GetChargeOffsets(i + 1);

                foreach (var offset in offsets)
                {
                    var offsetMz = mz + offset;
                    var peak = specMatch.Spectrum.FindPeak(offsetMz, _tolerance);
                    
                    if (peak != null)
                    {
                        peaks.Remove(peak);
                    }
                }
            }

            var spectrum = new Spectrum(peaks, specMatch.ScanNum);
            var filteredSpecMatch = new SpectrumMatch(specMatch.Peptide, spectrum, specMatch.ScanNum, charge);
            return filteredSpecMatch;
        }

        public List<SpectrumMatch> FilterMatches(List<SpectrumMatch> matches, bool filterNoise=true)
        {
            for (int i = 0; i < matches.Count; i++)
            {
                matches[i] = Filter(matches[i], filterNoise);
            }
            return matches;
        }
    }
}
