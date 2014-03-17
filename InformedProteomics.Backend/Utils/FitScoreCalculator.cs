using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices.ComTypes;

namespace InformedProteomics.Backend.Utils
{
    public class FitScoreCalculator
    {
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

        public static double GetFitNormalizedByTheoMaxIsotope(double[] theorPeakList, double[] observedPeakList)
        {
            if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 1.0;

            var maxIndex = -1;
            var maxIntensity = 0.0;
            for (var i = 0; i < theorPeakList.Length; i++)
            {
                if (theorPeakList[i] > maxIntensity)
                {
                    maxIndex = i;
                    maxIntensity = theorPeakList[i];
                }
            }

            var maxObs = observedPeakList[maxIndex];
            if (Math.Abs(maxObs) <= 0) return 1.0;

            var normalizedObs = new double[observedPeakList.Length];
            for (var i = 0; i < observedPeakList.Length; i++)
            {
                var normalizedValue = observedPeakList[i]/maxObs;
                normalizedObs[i] = normalizedValue > 1 ? 1 : normalizedValue;
            }

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

        // the larger the better
        public static double GetCosine(float[] theorPeakList, float[] observedPeakList)
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

            return innerProduct/Math.Sqrt(magnitudeTheo*magnitudeObs);
        }

        //public static double GetCosine2(float[] theorPeakList, float[] observedPeakList)
        //{
        //    if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 0;

        //    var innerProduct = theorPeakList.Select((t, i) => t * observedPeakList[i]).Sum();

        //    var magnitudeTheo = Math.Sqrt(theorPeakList.Sum(t => t * t));
        //    var magnitudeObs = Math.Sqrt(observedPeakList.Sum(t => t * t));

        //    return innerProduct / (magnitudeTheo * magnitudeObs);
        //}

        public static double GetFitOfNormalizedVectors(float[] normTheorPeakList, float[] normObservedPeakList)
        {
            if (normTheorPeakList.Length != normObservedPeakList.Length || normTheorPeakList.Length == 0) return 1.0;
            float sumSquareOfDiffs = 0;
            float sumSquareOfTheor = 0;
            for (var i = 0; i < normTheorPeakList.Length; i++)
            {
                var diff = normTheorPeakList[i] - normObservedPeakList[i];

                sumSquareOfDiffs += (diff * diff);
                sumSquareOfTheor += (normTheorPeakList[i] * normTheorPeakList[i]);
            }

            var fitScore = sumSquareOfDiffs / sumSquareOfTheor;
            if (float.IsNaN(fitScore) || fitScore > 1) fitScore = 1;

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

        static public float GetPearsonCorrelation(float[] v1, float[] v2)
        {
            var dimension = v1.Length;
            if (dimension <= 1 || dimension != v2.Length) return 0f;

            var m1 = 0f;
            var m2 = 0f;

            for (var i = 0; i < dimension; i++)
            {
                m1 += v1[i];
                m2 += v2[i];
            }

            m1 /= dimension;
            m2 /= dimension;

            var cov = 0f;
            var s1 = 0f;
            var s2 = 0f;

            for (var i = 0; i < dimension; i++)
            {
                var d1 = v1[i] - m1;
                var d2 = v2[i] - m2;
                cov += d1*d2;
                s1 += d1*d1;
                s2 += d2*d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0f : cov / (float)Math.Sqrt(s1 * s2);
        }
    }
}
