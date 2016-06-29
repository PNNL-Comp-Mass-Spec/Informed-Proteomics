using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics;

namespace InformedProteomics.Backend.Data.Composition
{
    public class IsoProfilePredictor
    {
        public const int MaxNumIsotopes = 100;
        public const double IsotopeRelativeIntensityThreshold = 0.1;

        public static IsotopomerEnvelope GetIsotopomerEnvelop(int c, int h, int n, int o, int s)
        {
            return Predictor.GetIsotopomerEnvelope(c, h, n, o, s);
        }

        public IsotopomerEnvelope GetIsotopomerEnvelope(int c, int h, int n, int o, int s)
        {
            var dist = new double[MaxNumIsotopes];
            var means = new double[_possibleIsotopeCombinations[0][0].Length + 1];
            var exps = new double[means.Length];
            for (var i = 0; i < means.Length; i++) // precalculate means and thier exps
            {
                means[i] = c * ProbC[i] + h * ProbH[i] + n * ProbN[i] + o * ProbO[i] + s * ProbS[i];
                exps[i] = Math.Exp(means[i]);
            }

            // This assumes that the envelop is unimodal.
            var maxHeight = 0.0;
            var isotopeIndex = 0;
            var mostIntenseIsotopomerIndex = -1;
            for (; isotopeIndex < MaxNumIsotopes; isotopeIndex++)
            {
                foreach (var isopeCombinations in _possibleIsotopeCombinations[isotopeIndex])
                {
                    dist[isotopeIndex] += GetIsotopeProbability(isopeCombinations, means, exps);
                }
                if (Double.IsInfinity(dist[isotopeIndex]))
                {
                    throw new NotFiniteNumberException();
                }
                if (dist[isotopeIndex] > maxHeight)
                {
                    maxHeight = dist[isotopeIndex];
                    mostIntenseIsotopomerIndex = isotopeIndex;
                }
                else if (dist[isotopeIndex] < maxHeight * IsotopeRelativeIntensityThreshold)
                {
                    break;
                }
            }

            var truncatedDist = new double[isotopeIndex];
            for (var i = 0; i < isotopeIndex; i++)
            {
                truncatedDist[i] = dist[i] / maxHeight;
            }

            return new IsotopomerEnvelope(truncatedDist, mostIntenseIsotopomerIndex);
        }

        public static IsoProfilePredictor Predictor;

        static IsoProfilePredictor()
        {
            Predictor = new IsoProfilePredictor();
        }

        public IsoProfilePredictor()
        {
            ProbC = DefaultProbC;
            ProbH = DefaultProbH;
            ProbN = DefaultProbN;
            ProbO = DefaultProbO;
            ProbS = DefaultProbS;
            ComputePossibleIsotopeCombinations(MaxNumIsotopes);
        }

        public IsoProfilePredictor(double[] probC, double[] probH, double[] probN, double[] probO, double[] probS)
        {
            ProbC = probC;
            ProbH = probH;
            ProbN = probN;
            ProbO = probO;
            ProbS = probS;
            ComputePossibleIsotopeCombinations(MaxNumIsotopes);
        }

        public static readonly double[] DefaultProbC = { .9893, 0.0107, 0, 0 };
        public static readonly double[] DefaultProbH = { .999885, .000115, 0, 0 };
        public static readonly double[] DefaultProbN = { 0.99632, 0.00368, 0, 0 };
        public static readonly double[] DefaultProbO = { 0.99757, 0.00038, 0.00205, 0 };
        public static readonly double[] DefaultProbS = { 0.9493, 0.0076, 0.0429, 0.0002 };


        public double[] ProbC { get; private set; }
        public double[] ProbH { get; private set; }
        public double[] ProbN { get; private set; }
        public double[] ProbO { get; private set; }
        public double[] ProbS { get; private set; }

        private int[][][] _possibleIsotopeCombinations;

        private void ComputePossibleIsotopeCombinations(int max) // called just once. 
        {
            var comb = new List<int[]>[max + 1];
            var maxIsotopeNumberInElement = ProbC.Length - 1;
            comb[0] = new List<int[]> { (new int[maxIsotopeNumberInElement]) };

            for (var n = 1; n <= max; n++)
            {
                comb[n] = new List<int[]>();
                for (var j = 1; j <= maxIsotopeNumberInElement; j++)
                {
                    var index = n - j;
                    if (index < 0) continue;
                    foreach (var t in comb[index])
                    {
                        var add = new int[maxIsotopeNumberInElement];
                        add[j - 1]++;
                        for (var k = 0; k < t.Length; k++)
                            add[k] += t[k];
                        var toAdd = comb[n].Select(v => !v.Where((t1, p) => t1 != add[p]).Any()).All(equal => !equal);
                        if (toAdd) comb[n].Add(add);
                    }
                }
            }
            var possibleIsotopeCombinations = new int[max][][];
            for (var i = 0; i < possibleIsotopeCombinations.Length; i++)
            {
                possibleIsotopeCombinations[i] = new int[comb[i].Count][];
                var j = 0;
                foreach (var t in comb[i])
                {
                    possibleIsotopeCombinations[i][j++] = t;
                }
            }
            _possibleIsotopeCombinations = possibleIsotopeCombinations;
        }

        private double GetIsotopeProbability(int[] number, double[] means, double[] exps)
        {
            var prob = 1.0;
            for (var i = 1; i <= Math.Min(ProbC.Length - 1, number.Length); i++)
            {
                var mean = means[i];
                var exp = exps[i];
                if (number[i - 1] == 0) prob *= exp;
                else
                    prob *=
                        (Math.Pow(mean, number[i - 1]) * exp / SpecialFunctions.Factorial(number[i - 1]));
            }
            return prob;
        }
    }
}
