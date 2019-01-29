using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.MathAndStats
{
    /// <summary>
    /// Static class with array handling utilities
    /// </summary>
    public static class ArrayUtil
    {
        /// <summary>
        /// Create a string to display the array values.
        /// </summary>
        /// <param name="array">The array</param>
        /// <param name="delimiter">Delimiter character</param>
        /// <param name="format">Optional. A string to use to format each value. Must contain the colon, so something like ':0.000'</param>
        public static string ToString<T>(T[] array, string delimiter = "\t", string format = "")
        {
            var s = new StringBuilder();
            var formatString = "{0" + format + "}";

            for (var i = 0; i < array.Length; i++)
            {
                if (i < array.Length - 1)
                {
                    s.AppendFormat(formatString + delimiter, array[i]);
                }
                else
                {
                    s.AppendFormat(formatString, array[i]);
                }
            }

            return s.ToString();
        }

        /// <summary>
        /// Create a string to display the array values.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array"></param>
        /// <param name="delimiter"></param>
        /// <param name="format"></param>
        /// <returns></returns>
        public static string ToString<T>(T[][] array, string delimiter = "\t", string format = "")
        {
            var s = new StringBuilder();
            var formatString = "{0" + format + "}";

            foreach (var item in array)
            {
                for (var i = 0; i < item.Length; i++)
                {
                    if (i < item.Length - 1)
                    {
                        s.AppendFormat(formatString + delimiter, item[i]);
                    }
                    else
                    {
                        s.AppendFormat(formatString, item[i]);
                    }
                }
                s.Append("\n");
            }

            return s.ToString();
        }

        /// <summary>
        /// Return an array ranking every value in <paramref name="values"/>, highest value 1 and lowest value n
        /// </summary>
        /// <param name="values"></param>
        /// <param name="lowerBoundValue"></param>
        /// <returns></returns>
        public static int[] GetRankings(IEnumerable<double> values, double lowerBoundValue = 0.0d)
        {
            var temp = new List<KeyValuePair<double, int>>();
            var i = 0;
            foreach (var v in values)
            {
                if (v > lowerBoundValue) temp.Add(new KeyValuePair<double, int>(v, i));
                i++;
            }

            var ranking = 1;
            var rankingList = new int[i];
            foreach (var t in temp.OrderByDescending(x => x.Key))
            {
                rankingList[t.Value] = ranking++;
            }
            return rankingList;
        }

        /// <summary>
        /// Return an array ranking every value in <paramref name="values"/>, highest value 1 and lowest value n
        /// </summary>
        /// <param name="values"></param>
        /// <param name="median"></param>
        /// <param name="lowerBoundValue"></param>
        /// <returns></returns>
        public static int[] GetRankings(IEnumerable<double> values, out double median, double lowerBoundValue = 0.0d)
        {
            var temp = new List<KeyValuePair<double, int>>();
            var i = 0;
            foreach (var v in values)
            {
                if (v > lowerBoundValue) temp.Add(new KeyValuePair<double, int>(v, i));
                i++;
            }

            var ranking = 1;
            var rankingList = new int[i];
            var medianRanking = (int)Math.Max(Math.Round(0.5*i), 1);
            median = 0;
            foreach (var t in temp.OrderByDescending(x => x.Key))
            {
                if (ranking == medianRanking) median = t.Key;
                rankingList[t.Value] = ranking++;
            }
            return rankingList;
        }

        /// <summary>
        /// Kadane's algorithm
        /// </summary>
        /// <param name="a"></param>
        /// <param name="start"></param>
        /// <param name="len"></param>
        /// <returns></returns>
        public static int MaxSumSubarray(IList<int> a, out int start, out int len)
        {
            start = 0;
            len = 1;
            var sum = a[0];

            var curStart = 0;
            var curLen = 1;
            var curSum = a[0];

            for (var i = 1; i < a.Count; i++)
            {
                if (a[i] >= curSum + a[i])
                {
                    curStart = i;
                    curLen = 1;
                    curSum = a[i];
                }
                else
                {
                    curLen++;
                    curSum += a[i];
                }

                if ((curSum <= sum) && (curSum != sum || curLen >= len) &&
                    (curSum != sum || curLen != len || curStart >= start)) continue;

                start = curStart;
                len = curLen;
                sum = curSum;
            }

            return sum;
        }

        /// <summary>
        /// Kadane's algorithm
        /// https://en.wikipedia.org/wiki/Maximum_subarray_problem
        /// </summary>
        /// <param name="a"></param>
        /// <param name="start"></param>
        /// <param name="len"></param>
        /// <returns></returns>
        public static double MaxSumSubarray(IList<double> a, out int start, out int len)
        {
            start = 0;
            len = 1;
            var sum = a[0];

            var curStart = 0;
            var curLen = 1;
            var curSum = a[0];

            for (var i = 1; i < a.Count; i++)
            {
                if (a[i] >= curSum + a[i])
                {
                    curStart = i;
                    curLen = 1;
                    curSum = a[i];
                }
                else
                {
                    curLen++;
                    curSum += a[i];
                }

                if ((curSum <= sum) && (Math.Abs(curSum - sum) > 1E-10 || curLen >= len) &&
                    (Math.Abs(curSum - sum) > 1E-10 || curLen != len || curStart >= start)) continue;

                start = curStart;
                len = curLen;
                sum = curSum;
            }

            return sum;
        }
    }
}
