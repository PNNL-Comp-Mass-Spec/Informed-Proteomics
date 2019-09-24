using System;
using System.Linq;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoreDistribution
    {
        public ScoreDistribution(int minScore, int maxScore)
        {
            MinScore = minScore;
            MaxScore = maxScore;
            _eValueDistribution = new double[MaxScore - MinScore];
            _scoreTotal = 0.0;
        }

        public int MinScore { get; }   // inclusive
        public int MaxScore { get; }   // exclusive

        public void SetEValue(int score, double eValue)
        {
            _eValueDistribution[score - MinScore] = eValue;
        }

        public void AddEValueDist(ScoreDistribution otherDistribution, int deltaScore, double weight)
        {
            if (otherDistribution == null) return;

            for (var index = Math.Max(otherDistribution.MinScore, MinScore - deltaScore); index < otherDistribution.MaxScore; index++)
            {
                var delEValue = otherDistribution._eValueDistribution[index - otherDistribution.MinScore]*weight;
                _eValueDistribution[index + deltaScore - MinScore] += delEValue;
                _scoreTotal += delEValue;
            }
        }

        public double GetSpectralEValue(double score)
        {
            return GetSpectralEValue((int)Math.Round(score));
        }

        public double GetSpectralEValue(int score)
        {
            var area = 0.0;

            var minIndex = (score >= MinScore) ? score - MinScore : 0;
            for (var index = minIndex; index < _eValueDistribution.Length; index++)
            {
                area += _eValueDistribution[index];
            }

            double spectralEValue;
            if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
            {
                if (area.Equals(0.0)) area = this._eValueDistribution.Last(v => v > 0);

                spectralEValue = area / _scoreTotal;
            }
            else
            {
                spectralEValue = area; // / _scoreTotal;
            }

            if (spectralEValue < double.Epsilon) return double.Epsilon; // to avoid underflow
            return spectralEValue;
        }

        private readonly double[] _eValueDistribution;

        private double _scoreTotal;
    }
}
