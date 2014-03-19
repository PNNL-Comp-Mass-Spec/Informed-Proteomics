using System;

using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class IonProbability: Probability
    {
        public IonType Ion { get; private set; }

        public IonProbability(IonType ion, int f=0, int t=0)
            : base(f, t)
        {
            Ion = ion;
        }


        public static IonProbability operator +(IonProbability l, IonProbability r)
        {
            var added = new IonProbability(l.Ion) {Found = l.Found + r.Found, Total = l.Total + r.Total};
            return added;
        }
    }
}
