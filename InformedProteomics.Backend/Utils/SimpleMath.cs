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
            return num / denom;
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

        static public double GetCorrelation(double[] v1, double[] v2)
        {
            if (v1.Length <= 1) return 0.0;

            /* var d = 0.0;
             var n = new double[2];

             for (var i = 0; i < v1.Length; i++)
             {
                 n[0] += v1[i]*v1[i];
                 n[1] += v2[i]*v2[i];
                 d += v1[i]*v2[i];
             }

             if (n[0] <= 0 || n[1] <= 0) return 0;
             return d/Math.Sqrt(n[0]*n[1]);
             //*/
            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = GetSampleVariance(v1, m1);
            var s2 = GetSampleVariance(v2, m2);
            if (s1 <= 0 || s2 <= 0) return 0;
            var div = Math.Sqrt(s1 * s2);
            var rho = v1.Select((t, i) => (float)((t - m1) * (v2[i] - m2) / div)).Sum();
            return Math.Min(Math.Max(0, rho / (v1.Length - 1)), 1);
            //*/
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
    }
}
