using System;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    internal class ScoreDistribution
    {
        public ScoreDistribution(int minScore, int maxScore)
        {
            MinScore = minScore;
            MaxScore = maxScore;
            _eValueDistribution = new double[MaxScore - MinScore];
        }

        public int MinScore { get; private set; }   // inclusive
        public int MaxScore { get; private set; }   // exclusive

        public void SetEValue(int score, double eValue)
        {
            _eValueDistribution[score - MinScore] = eValue;
        }

        public void AddEValueDist(ScoreDistribution otherDistribution, int deltaScore, double weight)
        {
            if (otherDistribution == null) return;
            for (var index = Math.Max(otherDistribution.MinScore, MinScore - deltaScore);
                index < otherDistribution.MaxScore;
                index++)
            {
                var delEValue = otherDistribution._eValueDistribution[index - otherDistribution.MinScore]*weight;
                _eValueDistribution[index + deltaScore - MinScore] += delEValue;
            }
        }

        public double GetSpectralEValue(int score)
        {
            var spectralEValue = 0.0;

            var minIndex = (score >= MinScore) ? score - MinScore : 0;
            for (var index = minIndex; index < _eValueDistribution.Length; index++)
            {
                spectralEValue += _eValueDistribution[index];
            }
            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        private readonly double[] _eValueDistribution;
    }
}
