using System;
using System.Linq;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Backend.Utils
{
    public class FitScoreCalculator
    {
        public static double GetBhattacharyyaDistance(double[] v1, double[] v2, int count = -1, int v1Index = 0, int v2Index = 0 )
        {
            if (count == -1) count = v1.Length;
            if (count == 0 || v1Index + count > v1.Length || v2Index + count > v2.Length) return 0.0;

            var s1 = 0d;
            var s2 = 0d;

            for (var i = 0; i < count; i++)
            {
                s1 += v1[i + v1Index];
                s2 += v2[i + v2Index];
            }

            if (!(s1 > 0) || !(s2 > 0)) return Double.PositiveInfinity;

            var bc = 0d;
            for (var i = 0; i < count; i++)
            {
                var p = v1[i + v1Index] / s1;
                var q = v2[i + v2Index] / s2;
                bc += Math.Sqrt(p * q);
            }

            return -Math.Log(bc);
        }

        public static double GetHyperGeometricPvalue(int n, int k, int n1, int k1, bool upperTailProb = true)
        {
            if (k >= n) return 1.0d;

            var pValue = Hypergeometric.CDF(n, k, n1, k1);
            if (upperTailProb) pValue = 1 - pValue;
            else pValue = Math.Min(pValue, 1 - pValue);

            return pValue;
        }

        public static double GetRankSumPvalue(double n, double n1, double r1, bool upperTailProb = true)
        {
            var n2 = n - n1;
            var u1 = n1 * n2 + n1 * (n1 + 1) * 0.5 - r1;

            var meanU = 0.5 * (n1 * n2);
            //var sigU = Math.Sqrt(n1*n2*(n1 + n2 + 1)/12);
            var logSigU = 0.5 * (Math.Log(n1) + Math.Log(n2) + Math.Log(n1 + n2 + 1) - Math.Log(12));
            var sigU = Math.Exp(logSigU);

            var pValue = Normal.CDF(meanU, sigU, u1);

            if (upperTailProb) pValue = 1 - pValue;
            else pValue = Math.Min(pValue, 1 - pValue);

            return Math.Abs(pValue); //negative tiny value
        }

        /*
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
        }*/

        public static Tuple<double, double> GetDistanceAndCorrelation(double[] v1, double[] v2, int count = -1, int v1Index = 0, int v2Index = 0)
        {
            if (count == -1) count = v1.Length;
            if (count == 0 || v1Index + count > v1.Length || v2Index + count > v2.Length) return new Tuple<double, double>(0d, 0d);

            var s1 = 0d;
            var s2 = 0d;

            for (var i = 0; i < count; i++)
            {
                s1 += v1[i + v1Index];
                s2 += v2[i + v2Index];
            }

            if (!(s1 > 0) || !(s2 > 0))
            {
                return new Tuple<double, double>(1d, 0d);
            }

            var m1 = s1 / count;
            var m2  = s2 / count;

            // compute Pearson correlation
            var cov = 0.0;
            var c1 = 0.0;
            var c2 = 0.0;
            var corr = 0d;
            var bc = 0d;

            for (var i = 0; i < count; i++)
            {
                var e1 = v1[v1Index + i];
                var e2 = v2[v2Index + i];
                var d1 = e1 - m1;
                var d2 = e2 - m2;

                cov += d1 * d2;
                c1 += d1 * d1;
                c2 += d2 * d2;

                var p = e1 / s1;
                var q = e2 / s2;
                bc += Math.Sqrt(p * q);
            }

            if (c1 > 0 && c2 > 0 && cov > 0)
            {
                corr = cov / Math.Sqrt(c1 * c2);
            }
            var dist = -Math.Log(bc);

            return new Tuple<double, double>(dist, corr);
        }


        public static double GetPearsonCorrelation(double[] v1, double[] v2, int count = -1, int v1Index = 0, int v2Index = 0)
        {
            if (count == -1) count = v1.Length;
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

        public static double GetDotProduct(double[] theorPeakList, double[] observedPeakList)
        {
            if (theorPeakList.Length != observedPeakList.Length || theorPeakList.Length == 0) return 0;

            var innerProduct = 0.0;
            for (var i = 0; i < theorPeakList.Length; i++)
            {
                var theo = theorPeakList[i];
                var obs = observedPeakList[i];
                innerProduct += theo * obs;
            }

            return innerProduct;
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
