using System;

namespace InformedProteomics.Backend.MathAndStats
{
    public class OtsuThreshold
    {
        // find otsu threshold
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

        private static double RunOtsuMethod(int[] hist, double minX, double intervalX)
        {
            var nBins = hist.Length;
            var vet = new double[nBins];
            double p1, p2, p12;
            for (var j = 1; j < nBins; j++)
            {
                //var th = minX + (j+1) * interval;
                p1 = Px(0, j, hist);
                p2 = Px(j + 1, nBins - 1, hist);

                p12 = p1 * p2;
                if (p12 < double.Epsilon) p12 = 1;

                var diff = (Mx(0, j, hist) * p2) - (Mx(j + 1, nBins - 1, hist) * p1);
                vet[j] = diff * diff / p12;
            }

            var thresholdIndex = findMax(vet);
            var threshold = minX + (thresholdIndex + 1) * intervalX;

            return threshold;
        }

        // function is used to compute the q values in the equation
        private static double Px(int init, int end, int[] hist)
        {
            var sum = 0d;
            int i;
            for (i = init; i <= end; i++)
                sum += hist[i];

            return sum;
        }

        // function is used to compute the mean values in the equation (mu)
        private static double Mx(int init, int end, int[] hist)
        {
            var sum = 0d;
            int i;
            for (i = init; i <= end; i++)
                sum += i * hist[i];

            return sum;
        }

        // finds the maximum element in a vector
        private static int findMax(double[] vec)
        {
            double maxVec = 0;
            int idx = 0;
            int i;

            for (i = 1; i < vec.Length; i++)
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
