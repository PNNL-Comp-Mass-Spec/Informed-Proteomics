using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.MathAndStats
{
    /// <summary>
    /// Class for calculating Otsu threshold
    /// </summary>
    public class OtsuThreshold
    {
        /// <summary>
        /// Calculate the Otsu threshold for the provided data
        /// </summary>
        /// <param name="data"></param>
        /// <param name="minX"></param>
        /// <param name="maxX"></param>
        /// <param name="intervalX"></param>
        /// <param name="xLb"></param>
        /// <param name="xUb"></param>
        /// <returns></returns>
        public static double GetThreshold(double[] data, double minX, double maxX, double intervalX, double xLb, double xUb)
        {
            var nBins = (int)Math.Ceiling((maxX - minX) / intervalX);
            var hist = new int[nBins];

            foreach (var x in data)
            {
                if (x < xLb || x > xUb) continue;

                var binIndex = (int)(Math.Ceiling((x - minX) / intervalX) - 1);
                if (binIndex < 0) binIndex = 0;
                else if (binIndex >= nBins) binIndex = nBins - 1;
                hist[binIndex]++;
            }

            return RunOtsuMethod(hist, minX, intervalX);
        }

        /// <summary>
        /// Calculate the Otsu threshold for the provided data
        /// </summary>
        /// <param name="data"></param>
        /// <param name="minX"></param>
        /// <param name="maxX"></param>
        /// <param name="intervalX"></param>
        /// <param name="minRow"></param>
        /// <param name="maxRow"></param>
        /// <param name="minCol"></param>
        /// <param name="maxCol"></param>
        /// <param name="xLb"></param>
        /// <param name="xUb"></param>
        /// <returns></returns>
        public static double GetThreshold(double[][] data, double minX, double maxX, double intervalX,
            int minRow, int maxRow, int minCol, int maxCol, double xLb, double xUb)
        {
            var nBins = (int)Math.Ceiling((maxX - minX) / intervalX);
            var hist = new int[nBins];

            for(var i = minRow; i < maxRow; i++)
            {
                for (var j = minCol; j < maxCol; j++)
                {
                    var x = data[i][j];
                    if (x < xLb || x > xUb) continue;

                    var binIndex = (int)(Math.Ceiling((x - minX) / intervalX) - 1);
                    if (binIndex < 0) binIndex = 0;
                    else if (binIndex >= nBins) binIndex = nBins - 1;

                    hist[binIndex]++;
                }
            }

            return RunOtsuMethod(hist, minX, intervalX);
        }

        private static double RunOtsuMethod(IReadOnlyList<int> hist, double minX, double intervalX)
        {
            var nBins = hist.Count;
            var vet = new double[nBins];
            for (var j = 1; j < nBins; j++)
            {
                //var th = minX + (j+1) * interval;
                var p1 = Px(0, j, hist);
                var p2 = Px(j + 1, nBins - 1, hist);

                var p12 = p1 * p2;
                if (p12 < double.Epsilon) p12 = 1;

                var diff = (Mx(0, j, hist) * p2) - (Mx(j + 1, nBins - 1, hist) * p1);
                vet[j] = diff * diff / p12;
            }

            var thresholdIndex = findMax(vet);
            var threshold = minX + (thresholdIndex + 1) * intervalX;

            return threshold;
        }

        // function is used to compute the q values in the equation
        private static double Px(int init, int end, IReadOnlyList<int> hist)
        {
            var sum = 0d;
            int i;
            for (i = init; i <= end; i++)
                sum += hist[i];

            return sum;
        }

        // function is used to compute the mean values in the equation (mu)
        private static double Mx(int init, int end, IReadOnlyList<int> hist)
        {
            var sum = 0d;
            int i;
            for (i = init; i <= end; i++)
                sum += i * hist[i];

            return sum;
        }

        // finds the maximum element in a vector
        private static int findMax(IReadOnlyList<double> vec)
        {
            double maxVec = 0;
            var idx = 0;
            int i;

            for (i = 1; i < vec.Count; i++)
            {
                if (vec[i] > maxVec)
                {
                    maxVec = vec[i];
                    idx = i;
                }
            }
            return idx;
        }
    }
}
