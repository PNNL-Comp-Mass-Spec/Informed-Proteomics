using System;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// IntRange: min, max, and test functions for working with ranges
    /// </summary>
    public class IntRange : IComparable<IntRange>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        public IntRange(int min, int max)
        {
            Min = min;
            Max = max;
        }

        /// <summary>
        /// Check if <paramref name="value"/> is within this range
        /// </summary>
        /// <param name="value"></param>
        /// <returns></returns>
        public bool Contains(int value)
        {
            return value >= Min && value <= Max;
        }

        /// <summary>
        /// Check if 2 IntRanges overlap
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Overlaps(IntRange other)
        {
            return Contains(other.Min) || Contains(other.Max);
        }

        /// <summary>
        /// Create a new IntRange that is the union of 2 other ranges
        /// </summary>
        /// <param name="range1"></param>
        /// <param name="range2"></param>
        /// <returns></returns>
        public static IntRange Union(IntRange range1, IntRange range2)
        {
            return new IntRange(Math.Min(range1.Min, range2.Min), Math.Max(range1.Max, range2.Max));
        }

        /// <summary>
        /// Add <paramref name="v"/> to the values covered by this IntRange
        /// </summary>
        /// <param name="v"></param>
        public void Add(int v)
        {
            if (v < Min)
            {
                Min = v;
            }

            if (v > Max)
            {
                Max = v;
            }
        }

        /// <summary>
        /// Minimum of this range
        /// </summary>
        public int Min { get; private set; }

        /// <summary>
        /// Maximum of this range
        /// </summary>
        public int Max { get; private set; }

        /// <summary>
        /// Length of this range
        /// </summary>
        public int Length => Max - Min + 1;

        /// <summary>
        /// Compare 2 IntRanges
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(IntRange other)
        {
            return Min.CompareTo(other.Min);
        }
    }
}