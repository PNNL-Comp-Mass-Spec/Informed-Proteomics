using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class FitScore
    {
        public double Intensity { get; private set; }
        public double Score { get; private set; }

        public FitScore(double intensity, double score)
        {
            Intensity = intensity;
            Score = score;
        }
    }

    public class CompareFitScoreByIntensity : IComparer<FitScore>
    {

        public int Compare(FitScore x, FitScore y)
        {
            return x.Intensity.CompareTo(y.Intensity);
        }
    }

    public class CompareFitScoreByScore : IComparer<FitScore>
    {
        public int Compare(FitScore x, FitScore y)
        {
            return x.Score.CompareTo(y.Score);
        }
    }
}
