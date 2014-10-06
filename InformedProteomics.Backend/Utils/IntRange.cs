using System;

namespace InformedProteomics.Backend.Utils
{
    public class IntRange : IComparable<IntRange>
    {
        public IntRange(int min, int max)
        {
            Min = min;
            Max = max;
        }

        public bool Contains(int value)
        {
            return value >= Min && value <= Max;
        }

        public bool Overlaps(IntRange other)
        {
            return Contains(other.Min) || Contains(other.Max);
        }

        public static IntRange Union(IntRange range1, IntRange range2)
        {
            return new IntRange(Math.Min(range1.Min, range2.Min), Math.Max(range1.Max, range2.Max));
        }

        public int Min { get; private set; }
        public int Max { get; private set; }

        public int CompareTo(IntRange other)
        {
            return Min.CompareTo(other.Min);
        }
    }
}