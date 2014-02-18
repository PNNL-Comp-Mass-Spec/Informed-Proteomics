using System;
using System.Linq;

namespace InformedProteomics.Backend.Utils
{
    public class FitScoreCalculator
    {
        public static double GetFit(double[] theorPeakList, double[] observedPeakList)
        {
            if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 1.0;

            var maxObs = observedPeakList.Max();
            if (Math.Abs(maxObs - 0) < float.Epsilon) maxObs = double.PositiveInfinity;
            var normalizedObs = observedPeakList.Select(p => p / maxObs).ToList();

            var maxTheor = theorPeakList.Max();
            var normalizedTheo = theorPeakList.Select(p => p / maxTheor).ToList();

            double sumSquareOfDiffs = 0;
            double sumSquareOfTheor = 0;
            for (var i = 0; i < normalizedTheo.Count; i++)
            {
                var diff = normalizedObs[i] - normalizedTheo[i];

                sumSquareOfDiffs += (diff * diff);
                sumSquareOfTheor += (normalizedTheo[i] * normalizedTheo[i]);
            }

            var fitScore = sumSquareOfDiffs / sumSquareOfTheor;
            if (double.IsNaN(fitScore) || fitScore > 1) fitScore = 1;

            return fitScore;
        }

        public static double GetFitOfNormalizedVectors(float[] normTheorPeakList, float[] normObservedPeakList)
        {
            if (normTheorPeakList.Length != normObservedPeakList.Length || normTheorPeakList.Length == 0) return 1.0;
            double sumSquareOfDiffs = 0;
            double sumSquareOfTheor = 0;
            for (var i = 0; i < normTheorPeakList.Length; i++)
            {
                var diff = normTheorPeakList[i] - normObservedPeakList[i];

                sumSquareOfDiffs += (diff * diff);
                sumSquareOfTheor += (normTheorPeakList[i] * normTheorPeakList[i]);
            }

            var fitScore = sumSquareOfDiffs / sumSquareOfTheor;
            if (double.IsNaN(fitScore) || fitScore > 1) fitScore = 1;

            return fitScore;
        }

    }
}
