using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Scoring.TopDown
{
    public class IcScores
    {
        public IcScores(int nMatchedFragments, double score, string modifications)
        {
            NumMatchedFrags = nMatchedFragments;
            Score = score;
            Modifications = modifications;
        }

        public int NumMatchedFrags { get; private set; }
        public double Score { get; private set; } // this score is used to calculate p-value by generating function

        public string Modifications { get; private set; }

        public override string ToString()
        {
            return string.Join("\t",
                new[]
                {
                    NumMatchedFrags, Score,
                });
        }

        public static string GetScoreNames()
        {
            return "#MatchedFragments\tScore";
        }
    }
}
