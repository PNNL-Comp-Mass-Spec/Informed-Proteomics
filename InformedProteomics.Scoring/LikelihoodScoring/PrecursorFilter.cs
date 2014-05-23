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

        public PrecursorFilter(Dictionary<int, PrecursorOffsets> offsets, Tolerance tolerance)
        {
            _offsets = offsets;
            _tolerance = tolerance;

            var maxCharge = _offsets.Keys.Max();

            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, maxCharge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        public PrecursorFilter(int maxCharge, Tolerance tolerance)
        {
            _offsets = new Dictionary<int, PrecursorOffsets>();
            _tolerance = tolerance;

            for (int i = 0; i < maxCharge; i++)
            {
                _offsets.Add(i+1, new PrecursorOffsets(i+1));
            }

            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, maxCharge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        public void SetChargeOffsets(PrecursorOffsets precursorOffsets)
        {
            var charge = precursorOffsets.Charge;
            _offsets[charge] = precursorOffsets;
        }

        private SpectrumMatch Filter(SpectrumMatch specMatch)
        {
            var charge = specMatch.PrecursorCharge;
            var peaks = specMatch.Spectrum.Peaks;
            var offsetLists = _offsets[charge];
            var indexes = new List<int>();

            for (int i = 0; i < offsetLists.Charge; i++)
            {
                var ion = _precursorIonTypes[i].GetIon(specMatch.PeptideComposition);
                var mz = ion.GetMonoIsotopicMz();

                var offsets = offsetLists.GetChargeOffsets(i + 1);
/*                for (int j=0; j < peaks.Length; j++)
                {
                    var peak = peaks[j];
                    if (peak.Mz >= (mz - 100/(i+1)) && peak.Mz <= (mz + 100/(i+1)))
                    {
                        indexes.Add(j);
                    }
                } */
                
                foreach (var offset in offsets)
                {
                    var offsetMz = mz + offset;
                    var peakIndex = specMatch.Spectrum.FindPeakIndex(offsetMz, _tolerance);
                    if (peakIndex > 0) indexes.Add(peakIndex);
                }
            }
            indexes = indexes.Distinct().ToList();
            var filteredPeaks = peaks.Where((t, i) => !indexes.Contains(i)).ToList();

            var spectrum = new Spectrum(filteredPeaks, specMatch.ScanNum);
            var filteredSpecMatch = new SpectrumMatch(specMatch.Peptide, spectrum, specMatch.ScanNum, charge, specMatch.Sequence);
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
    }
}
