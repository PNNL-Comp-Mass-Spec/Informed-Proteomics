using System;
using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class Probability <T1>
    {
        public T1 DataLabel { get; private set; }
        public int Found { get; set; }
        public int Total { get; set; }

        public double Prob
        {
            get { return Math.Round((double)(Found) / (Total), 5); }
        }
        public Probability(T1 dataLabel, int f=0, int t=0)
        {
            DataLabel = dataLabel;
            Found = f;
            Total = t;
        }


        public static Probability<T1> operator +(Probability<T1> l, Probability<T1> r)
        {
            var added = new Probability<T1>(l.DataLabel, l.Found + r.Found, l.Total + r.Total);
            return added;
        }
    }

    public class CompareByProbability<T1> : IComparer<Probability<T1>>
    {
        public int Compare(Probability<T1> x, Probability<T1> y)
        {
            return x.Prob.CompareTo(y.Prob);
        }
    }
}
