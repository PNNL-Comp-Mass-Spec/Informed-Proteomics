using System;
using System.Linq;

namespace InformedProteomics.Backend.Utils
{
    public class FitScoreCalculator
    {
        public static double GetPearsonCorrelation(double[] v1, double[] v2)
        {
            var dimension = v1.Length;
            if (dimension == 0 || dimension != v2.Length) return 0.0;
            if (dimension == 1) return 1.0;

            // Compute means
            var m1 = 0.0;
            var m2 = 0.0;

            for (var i = 0; i < dimension; i++)
            {
                m1 += v1[i];
                m2 += v2[i];
            }

            m1 /= dimension;
            m2 /= dimension;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < dimension; i++)
            {
                var d1 = v1[i] - m1;
                var d2 = v2[i] - m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0f : cov / Math.Sqrt(s1 * s2);
        }

        public static double GetPearsonCorrelation(double[] v1, int v1Index, double[] v2, int v2Index, int count)
        {
            if (count == 0 || v1Index + count > v1.Length || v2Index + count > v2.Length) return 0.0;
            if (count == 1) return 1.0;

            // Compute means
            var m1 = 0.0;
            var m2 = 0.0;

            for (var i = 0; i < count; i++)
            {
                m1 += v1[v1Index + i];
                m2 += v2[v2Index + i];
            }

            m1 /= count;
            m2 /= count;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < count; i++)
            {
                var d1 = v1[v1Index + i] - m1;
                var d2 = v2[v2Index + i] - m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0f : cov / Math.Sqrt(s1 * s2);
        }


        // the larger the better
        public static double GetCosine(double[] theorPeakList, double[] observedPeakList)
        {
            if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 0;

            var innerProduct = 0.0;
            var magnitudeTheo = 0.0;
            var magnitudeObs = 0.0;
            for (var i = 0; i < theorPeakList.Length; i++)
            {
                var theo = theorPeakList[i];
                var obs = observedPeakList[i];
                innerProduct += theo * obs;
                magnitudeTheo += theo * theo;
                magnitudeObs += obs * obs;
            }

            return innerProduct / Math.Sqrt(magnitudeTheo * magnitudeObs);
        }

        // the smaller the better
        public static double GetDeconToolsFit(double[] theorPeakList, double[] observedPeakList)
        {
            if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 1.0;

            var maxObs = observedPeakList.Max();
            if (Math.Abs(maxObs - 0) < float.Epsilon) maxObs = double.PositiveInfinity;
            var normalizedObs = observedPeakList.Select(p => p / maxObs).ToList();

            double sumSquareOfDiffs = 0;
            double sumSquareOfTheor = 0;
            for (var i = 0; i < theorPeakList.Length; i++)
            {
                var diff = normalizedObs[i] - theorPeakList[i];

                sumSquareOfDiffs += (diff * diff);
                sumSquareOfTheor += (theorPeakList[i] * theorPeakList[i]);
            }

            var fitScore = sumSquareOfDiffs / sumSquareOfTheor;
            if (double.IsNaN(fitScore) || fitScore > 1) fitScore = 1;

            return fitScore;
        }

        public static double GetFitOfNormalizedVectors(double[] normTheorPeakList, double[] normObservedPeakList)
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
