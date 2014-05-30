using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class PrecursorFilter
    {
        public PrecursorFilter(PrecursorOffsets offsets, Tolerance tolerance)
        {
            _offsets = offsets;
            _tolerance = tolerance;

            var maxCharge = offsets.Charge;

            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, maxCharge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        public SpectrumMatch Filter(SpectrumMatch specMatch)
        {
            var charge = specMatch.PrecursorCharge;
            var peaks = specMatch.Spectrum.Peaks;
            var indexes = new List<int>();

            for (int i = 0; i < _offsets.Charge; i++)
            {
                var ion = _precursorIonTypes[i].GetIon(specMatch.PeptideComposition);
                var mz = ion.GetMonoIsotopicMz();

                var offsets = _offsets.GetChargeOffsets(i + 1);
                
                foreach (var offset in offsets)
                {
                    var offsetMz = mz + offset;
                    var peakIndex = specMatch.Spectrum.FindPeakIndex(offsetMz, _tolerance);
                    if (peakIndex > 0) indexes.Add(peakIndex);
                }
            }
            indexes = indexes.Distinct().ToList();
            var filteredPeaks = peaks.Where((t, i) => !indexes.Contains(i)).ToList();

            var spectrum = new Spectrum(filteredPeaks, specMatch.Spectrum.ScanNum);
            var filteredSpecMatch = new SpectrumMatch(specMatch.Sequence, spectrum, charge, specMatch.Decoy);
            return filteredSpecMatch;
        }

        public List<SpectrumMatch> FilterMatches(List<SpectrumMatch> matches)
        {
            for (int i = 0; i < matches.Count; i++)
            {
                matches[i] = Filter(matches[i]);
            }
            return matches;
        }
        private readonly PrecursorOffsets _offsets;
        private readonly List<IonType> _precursorIonTypes;
        private readonly Tolerance _tolerance;
    }
}
