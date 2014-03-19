using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class FitScoreList: List<FitScore>
    {
        public FitScoreList(double[] intensities, double[] scores)
        {
            int counter = 0;
            if (intensities != null && intensities.Length > 0)
                counter = intensities.Length;
            else if (scores != null && scores.Length > 0)
                counter = scores.Length;

            for (int i = 0; i < counter; i++)
            {
                double intensity = 0.0;
                double score = 0.0;

                if (intensities != null)
                    intensity = intensities[i];
                if (scores != null)
                    score = scores[i];

                Add(new FitScore(intensity, score));
            }
        }

        public FitScoreList(IEnumerable<FitScore> fitscores)
        {
            foreach (var fitscore in fitscores)
            {
                Add(fitscore);
            }
        }

        public FitScoreList()
        {
        }

        public double[] Scores
        {
            get
            {
                var scores = new double[Count];
                for (int i = 0; i < Count; i++)
                {
                    scores[i] = this[i].Score;
                }
                return scores;
            }
        }

        public double[] Intensities
        {
            get
            {
                var intensities = new double[Count];
                for (int i = 0; i < Count; i++)
                {
                    intensities[i] = this[i].Intensity;
                }
                return intensities;
            }
        }

        public FitScore MaxScore(ScoreMethod method)
        {
            FitScore score = null;
            double bestScore = 0.0;

            foreach (var sc in this)
            {
                double scScore = sc.Score;
                if (method == ScoreMethod.FitScore)
                    scScore = 1 - scScore;

                if (scScore >= bestScore)
                {
                    score = sc;
                    bestScore = scScore;
                }
            }
            return score;
        }
    }
}
