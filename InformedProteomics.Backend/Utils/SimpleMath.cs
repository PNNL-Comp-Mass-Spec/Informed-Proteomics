using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Utils
{
    public class SimpleMath
    {
        public static int NChooseK(int n, int k)
        {
            int num = 1;
            for (int i = 0; i < k; i++)
                num *= n - i;
            int denom = 1;
            for (int i = 0; i < k; i++)
                denom *= k - i;
            return num/denom;
        }

        private static readonly Dictionary<Tuple<int, int>, double> LogCombinations = new Dictionary<Tuple<int, int>, double>();

        public static double GetCombination(int n, int k)
        {
            var sum = GetLogCombination(n, k);
            return Math.Exp(sum);
        }

        public static double GetLogCombination(int n, int k)
        {
            if (LogCombinations.TryGetValue(new Tuple<int, int>(n, k), out var sum)) return sum;

            for (var i = 0; i < k; i++)
            {
                sum += Math.Log(n - i);
                sum -= Math.Log(i + 1);
            }

            LogCombinations.Add(new Tuple<int, int>(n, k), sum);
            return sum;
        }

        public static int[][] GetCombinationsWithRepetition(int n, int length)
        {
            if (n <= 0)
                return null;

            if (length == 0)
            {
                return new[] { new int[0] };
            }
            if (length == 1)
            {
                var combinations = new int[n][];
                for (int i = 0; i < n; i++)
                {
                    combinations[i] = new[] { i };
                }
                return combinations;
            }
            else
            {
                var prevCombinations = GetCombinationsWithRepetition(n, length - 1);
                var combinations = new List<int[]>();
                foreach (var combination in prevCombinations)
                {
                    int lastValue = combination.Last();
                    for (int j = lastValue; j < n; j++)
                    {
                        var newCombination = new int[combination.Length + 1];
                        Array.Copy(combination, newCombination, combination.Length);
                        newCombination[newCombination.Length - 1] = j;
                        combinations.Add(newCombination);
                    }
                }
                return combinations.ToArray();
            }
        }

        public static int[][] GetNtoTheKCombinations(int n, int length)
        {
            if (n <= 0)
                return null;

            if (length == 0)
            {
                return new[] { new int[0] };
            }
            if (length == 1)
            {
                var combinations = new int[n][];
                for (var i = 0; i < n; i++)
                {
                    combinations[i] = new[] { i };
                }
                return combinations;
            }
            else
            {
                var prevCombinations = GetNtoTheKCombinations(n, length - 1);
                var combinations = new List<int[]>();
                foreach (var combination in prevCombinations)
                {
                    for (var j = 0; j < n; j++)
                    {
                        var newCombination = new int[combination.Length + 1];
                        Array.Copy(combination, newCombination, combination.Length);
                        newCombination[newCombination.Length - 1] = j;
                        combinations.Add(newCombination);
                    }
                }
                return combinations.ToArray();
            }
        }

        static public double GetCorrelation(double[] v1, double[] v2)
        {
            if (v1.Length <= 1) return 0.0;

            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = GetSampleVariance(v1, m1);
            var s2 = GetSampleVariance(v2, m2);
            if (s1 <= 0 || s2 <= 0) return 0;
            var div = Math.Sqrt(s1 * s2);
            var rho = v1.Select((t, i) => (float)((t - m1) * (v2[i] - m2) / div)).Sum();
            return Math.Min(Math.Max(0, rho / (v1.Length - 1)), 1);
        }

        static public double GetSampleMean(double[] x)
        {
            var m = x.Sum();
            return m / x.Length;
        }

        static public double GetSampleVariance(double[] x, double m)
        {
            var var = x.Sum(v => (v - m) * (v - m));
            return var / (x.Length - 1);
        }

        static public float GetSampleMean(float[] x)
        {
            return x.Average();
        }

        static public float GetSampleVariance(float[] x, float m)
        {
            var variance = 0f;
            var length = x.Length;
            for (var i = 0; i < length; i++)
            {
                var diff = x[i] - m;
                variance += diff*diff;
            }
            return variance / (length - 1);
        }

        static public IEnumerable<int> GetFactors(int x)
        {
            for (var factor = 1; factor * factor <= x; factor++)
            {
                if (x % factor != 0) continue;
                yield return factor;
                if (factor * factor != x) yield return x / factor;
            }
        }

        static public double GetKLDivergence(double[] P, double[] Q)
        {
            var ret = P.Select((t, i) => t*(Math.Log(t) - Math.Log(Q[i]))).Sum();

            return ret;
        }
    }
}
