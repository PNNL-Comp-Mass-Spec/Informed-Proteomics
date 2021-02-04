using System;

namespace InformedProteomics.Scoring.BottomUp
{
    public class Partition : IComparable<Partition>
    {
        protected bool Equals(Partition other)
        {
            return Charge == other.Charge && NeutralPeptideMass.Equals(other.NeutralPeptideMass) && SegmentIndex == other.SegmentIndex;
        }

        public Partition(int charge, float neutralPeptideMass, int segmentIndex)
        {
            Charge = charge;
            NeutralPeptideMass = neutralPeptideMass;
            SegmentIndex = segmentIndex;
        }

        public int Charge { get; }
        public float NeutralPeptideMass { get; }
        public int SegmentIndex { get; }

        public int CompareTo(Partition other)
        {
            if (Charge < other.Charge)
            {
                return -1;
            }

            if (Charge > other.Charge)
            {
                return 1;
            }

            if (SegmentIndex < other.SegmentIndex)
            {
                return -1;
            }

            if (SegmentIndex > other.SegmentIndex)
            {
                return 1;
            }

            if (NeutralPeptideMass < other.NeutralPeptideMass)
            {
                return -1;
            }

            return NeutralPeptideMass > other.NeutralPeptideMass ? 1 : 0;
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
            {
                return false;
            }

            if (ReferenceEquals(this, obj))
            {
                return true;
            }

            if (obj.GetType() != GetType())
            {
                return false;
            }

            return Equals((Partition)obj);
        }

        public override int GetHashCode()
        {
            var hashCode = Charge;
            hashCode = (hashCode * 397) ^ NeutralPeptideMass.GetHashCode();
            hashCode = (hashCode * 397) ^ SegmentIndex;
            return hashCode;
        }
    }
}
