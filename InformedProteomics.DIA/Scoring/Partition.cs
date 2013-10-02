using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.DIA.Scoring
{
    public class Partition: IComparable<Partition>
    {
        public Partition(int charge, float neutralPeptideMass, int segmentIndex)
        {
            Charge = charge;
            NeutralPeptideMass = neutralPeptideMass;
            SegmentIndex = segmentIndex;
        }

        public int Charge { get; private set; }
        public float NeutralPeptideMass { get; private set; }
        public int SegmentIndex { get; private set; }


        public int CompareTo(Partition other)
        {
            if (Charge < other.Charge)
                return -1;
            if (Charge > other.Charge)
                return 1;
            if (SegmentIndex < other.SegmentIndex)
                return -1;
            if (SegmentIndex > other.SegmentIndex)
                return 1;
            if (NeutralPeptideMass < other.NeutralPeptideMass)
                return -1;
            return NeutralPeptideMass > other.NeutralPeptideMass ? 1 : 0;
        }
    }
}
