using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class CleavageOffsetFrequencyTable: OffsetFrequencyTable
    {
        public CleavageOffsetFrequencyTable(int searchWidth, int charge = 1, double binWidth = 1.005):
            base(searchWidth, charge, binWidth)
        {
        }

        public override void AddMatches(List<SpectrumMatch> matches)
        {
            throw new NotImplementedException();
        }
    }
}
