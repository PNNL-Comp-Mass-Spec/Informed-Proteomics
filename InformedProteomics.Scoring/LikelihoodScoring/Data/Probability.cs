﻿using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class Probability <T1>
    {
        public T1 Label { get; private set; }
        public double Found { get; set; }
        public double Total { get; set; }

        public double Prob
        {
            get { return (Found / Total); }
        }
        public Probability(T1 label, double found=0, double total=0)
        {
            Label = label;
            Found = found;
            Total = total;
        }


        public static Probability<T1> operator +(Probability<T1> l, Probability<T1> r)
        {
            var added = new Probability<T1>(l.Label, l.Found + r.Found, l.Total + r.Total);
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
