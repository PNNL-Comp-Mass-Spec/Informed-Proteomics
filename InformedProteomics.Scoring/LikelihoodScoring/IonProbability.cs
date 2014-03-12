using System;

using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class IonProbability
    {
        public int Found { get; set; }
        public int Total { get; set; }
        public IonType Ion { get; private set; }

        public double Probability
        {
            get { return Math.Round((double)(Found) / (Total), 5); }
        }
        public IonProbability(IonType ion, int f=0, int t=0)
        {
            Found = f;
            Total = t;
            Ion = ion;
        }


        public static IonProbability operator +(IonProbability l, IonProbability r)
        {
            var added = new IonProbability(l.Ion) {Found = l.Found + r.Found, Total = l.Total + r.Total};
            return added;
        }
    }
}
