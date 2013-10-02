using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.DIA.Scoring
{
    public class PrecursorOffsetFrequency
    {
        public PrecursorOffsetFrequency(int reducedCharge, float offset, float frequency, Tolerance tolerance)
        {
            ReducedCharge = reducedCharge;
            Offset = offset;
            Frequency = frequency;
            Tolerance = tolerance;
        }

        public int ReducedCharge { get; private set; }
        public float Offset { get; private set; }
        public float Frequency { get; private set; }
        public Tolerance Tolerance { get; private set; }
    }
}
