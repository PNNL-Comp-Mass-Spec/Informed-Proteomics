using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class FitScore
    {
        public double Intensity { get; }
        public double Score { get; }

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
            if (x == null)
            {
                return y == null ? 0 : -1;
            }

            if (y == null)
            {
                return 1;
            }

            return x.Intensity.CompareTo(y.Intensity);
        }
    }

    public class CompareFitScoreByScore : IComparer<FitScore>
    {
        public int Compare(FitScore x, FitScore y)
        {
            if (x == null)
            {
                return y == null ? 0 : -1;
            }

            if (y == null)
            {
                return 1;
            }

            return x.Score.CompareTo(y.Score);
        }
    }
}
