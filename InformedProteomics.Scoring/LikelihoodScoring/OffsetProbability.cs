using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class OffsetProbability: Probability
    {
        public double Offset { get; private set; }

        public OffsetProbability(double offset, int found = 0, int total = 0)
            : base(found, total)
        {
            Offset = offset;
        }


        public static OffsetProbability operator +(OffsetProbability l, OffsetProbability r)
        {
            var added = new OffsetProbability(l.Offset) { Found = l.Found + r.Found, Total = l.Total + r.Total };
            return added;
        }
    }

    public class CompareOffsteProbabilityByOffset : IComparer<OffsetProbability>
    {

        public int Compare(OffsetProbability x, OffsetProbability y)
        {
            return x.Offset.CompareTo(y.Offset);
        }
    }
}
