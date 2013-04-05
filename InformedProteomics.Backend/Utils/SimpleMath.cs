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
                int[][] prevCombinations = GetCombinationsWithRepetition(n, length - 1);
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
    }
}
