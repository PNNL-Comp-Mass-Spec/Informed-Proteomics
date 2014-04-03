using System;
using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class Probability
    {
        public int Found { get; set; }
        public int Total { get; set; }

        public double Prob
        {
            get { return Math.Round((double)(Found) / (Total), 5); }
        }
        public Probability(int f=0, int t=0)
        {
            Found = f;
            Total = t;
        }


        public static Probability operator +(Probability l, Probability r)
        {
            var added = new Probability(l.Found + r.Found, l.Total + r.Total);
            return added;
        }
    }

    public class CompareByProbability : IComparer<Probability>
    {
        public int Compare(Probability x, Probability y)
        {
            return x.Prob.CompareTo(y.Prob);
        }
    }
}
