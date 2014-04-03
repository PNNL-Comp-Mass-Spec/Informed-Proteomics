using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class PrecursorOffsetFrequencyTable: OffsetFrequencyTable
    {
        private readonly List<IonType> _precursorIonTypes;

        public PrecursorOffsetFrequencyTable(int searchWidth, int charge = 1, double binWidth = 1.005):
            base(searchWidth, charge, binWidth)
        {
            BaseIonType[] baseIons = { BaseIonType.Y };
            NeutralLoss[] neutralLosses = { NeutralLoss.NoLoss };

            var ionTypeFactory = new IonTypeFactory(baseIons, neutralLosses, charge);

            _precursorIonTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
        }

        public override void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                Total++;
                var ion = _precursorIonTypes[Charge - 1].GetIon(match.PeptideComposition);
                var monoIsotopicMz = ion.GetMonoIsotopicMz();
                var min = GetMinMz(monoIsotopicMz);
                var max = GetMaxMz(monoIsotopicMz);

                var peaks = match.Spectrum.Peaks;
                var offsetMzCollection = new List<double>();

                foreach (var peak in peaks)
                {
                    if (peak.Mz >= min && peak.Mz <= max)
                    {
                        offsetMzCollection.Add(peak.Mz - monoIsotopicMz);
                    }
                }
                AddOffsets(offsetMzCollection);
            }
        }
    }
}
