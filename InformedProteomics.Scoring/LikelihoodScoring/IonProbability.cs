using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class IonProbability: Probability
    {
        public string IonName { get; private set; }

        public IonProbability(string ion, int f=0, int t=0)
            : base(f, t)
        {
            IonName = ion;
        }


        public static IonProbability operator +(IonProbability l, IonProbability r)
        {
            var added = new IonProbability(l.IonName) {Found = l.Found + r.Found, Total = l.Total + r.Total};
            return added;
        }
    }
}
