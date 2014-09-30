using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Backend.SequenceTag
{
    public class ScoringFunction
    {
        public int[] PeakRankings { get; private set; }
        public List<DeconvolutedPeak> DeconvolutedPeaks;
        private static readonly Normal Gaussian = new Normal();

        public ScoringFunction(List<DeconvolutedPeak> _deconvolutedPeaks)
        {
            DeconvolutedPeaks = _deconvolutedPeaks;
            // compute peak rankings for scoring
            var sorted = DeconvolutedPeaks.Select((x, i) => new KeyValuePair<DeconvolutedPeak, int>(x, i)).OrderByDescending(x => x.Key.Intensity).ToList();
            PeakRankings = new int[sorted.Count];
            for (var ranking = 1; ranking < sorted.Count; ranking++)
            {
                PeakRankings[sorted[ranking - 1].Value] = ranking;
            }            
        }

        public double ComputeScoreByRankSumTest(List<int> matchedPeakIndex)
        {
            var n = (double)matchedPeakIndex.Count;
            var N = (double)PeakRankings.Length;
            double rankSum = 0;
            foreach (var i in matchedPeakIndex) rankSum += PeakRankings[i];

            double U = n * (N - n) + n * (n + 1) / 2 - rankSum;
            double zScore = (U - n * (N - n) / 2) / Math.Sqrt(n * (N - n) * (N + 1) / 12);

            return -Math.Log(1 - Gaussian.CumulativeDistribution(zScore), 10);
        }

        public double ComputeScoreByHyperGeometricTest(List<int> matchedPeakIndex)
        {
            double minMass = DeconvolutedPeaks[matchedPeakIndex[0]].Mass;
            double maxMass = DeconvolutedPeaks[matchedPeakIndex[matchedPeakIndex.Count - 1]].Mass;
            double massRange = maxMass - minMass;

            const double MassTolerance = 0.01; // need to fix this!!! dynamically assign depending on the m/z range of interest

            var n = (int)Math.Round(massRange / MassTolerance); // # of possible ion locations
            var k = matchedPeakIndex.Count; // # of theretical ions

            var n1 = matchedPeakIndex[matchedPeakIndex.Count - 1] - matchedPeakIndex[0] + 1; //# of observed ions
            var m = matchedPeakIndex.Count; // # of matched ions

            double score = NegativeLog10_NChooseK(n1, m) + NegativeLog10_NChooseK(n - n1, k - m) - NegativeLog10_NChooseK(n, k);

            return score;
        }

        public static double NegativeLog10_NChooseK(int n, int k)
        {
            double num = 0;
            for (double i = 0; i < k; i++)
                num -= Math.Log(n - i);

            double denom = 0;
            for (int i = 0; i < k; i++)
                denom -= Math.Log(k - i, 10);

            return num - denom;
        }

    }
}
