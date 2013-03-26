using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Single;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Scoring
{
    class FragmentXICSetScorer
    {
        public FragmentXICSet FragmentXICSet { get; private set; }
        public FragmentParameter Parameter { get; private set; }
        public float RawScore { get; private set; }
        public float Score { get; private set; }
        private readonly int _apexIndex;
        private readonly List<IonType> _usedIonTypes;
        private readonly IonType _precursorIon;

        public FragmentXICSetScorer(FragmentXICSet fragmentXICSet, int apexIndex, FragmentParameter par)
        {
            FragmentXICSet = fragmentXICSet;
            _apexIndex = apexIndex;
            Parameter = par;
            _usedIonTypes = new List<IonType>();
            foreach (var ion in fragmentXICSet.Keys)
            {
                if (ion.IsPrecursor) _precursorIon = ion;
                else _usedIonTypes.Add(ion);
            }
            Score = GetScore();
        }

        private float GetScore()
        {
            if (_usedIonTypes.Count == 0) return 0;
            var r = GetCorrelationMatrices(_usedIonTypes, Parameter);
            var xy = GetIonXICMatrices(FragmentXICSet, _usedIonTypes, _precursorIon, _apexIndex);
            RawScore = new MultipleCorrelationCoefficient(xy[0], xy[1], r[0], r[1]).Get();
            return SubScoreFactory.GetProductIonXICLikelihoodRatioScore(RawScore, _usedIonTypes.Count, Parameter);
        }

        private static List<DenseMatrix> GetIonXICMatrices(FragmentXICSet set, IList<IonType> ions, IonType precursorIon, int apexIndex)
        {
            var x = new DenseMatrix(set[ions[0]].Length, ions.Count);
            var y = new DenseMatrix(set[ions[0]].Length, ions.Count);

            var d = set.GetIonXIC(precursorIon);
            var nd = Divide(d, d[apexIndex]);
            
            var m = GetSampleMean(nd);
            var v = GetSampleVariance(nd, m);
            var nds = Standardize(nd, m, v);
            for (var j = 0; j < nds.Length; j++)
                y.At(j, 0, nds[j]);

            for (var i=0;i<ions.Count;i++)
            {
                var c = set.GetIonXIC(ions[i]);
                var nc = Divide(c, c[apexIndex]);
                var ncs = Standardize(nc, m, v);

                for (var j = 0; j < ncs.Length; j++)
                    x.At(j, i, ncs[j]);
            }

            var ret = new List<DenseMatrix> {x, y};
            return ret;
        }

        private static float GetSampleMean(float[] x)
        {
            var m = x.Sum();
            return m / x.Length;
        }

        private static float GetSampleVariance(float[] x, float m)
        {
            var var = x.Sum(v => (v - m) * (v - m));
            return var / (x.Length - 1);
        }


        private static DenseMatrix[] GetCorrelationMatrices(List<IonType> ions, FragmentParameter par)
        {
            var j = ions.Count;
            var rxx = new DenseMatrix(j, j, 0);
            var rxy = new DenseMatrix(j, 1, 0);
 
            for (var i = 0; i < j; i++)
            {
                for (var k = 0; k < j; k++)
                {
                    rxx.At(i, k, SubScoreFactory.GetIonXICCorrelationCoefficient(ions[i], ions[k], par));
                }
                rxy.At(i, 0, SubScoreFactory.GetIonXICCorrelationCoefficient(ions[i], par));
            }

            return new[] { rxx, rxy };
        }

        private static float[] Standardize(float[] x, float mean, float variance) // return standardized x (i.e., mean = 0, var = 1) 
        {
            var sx = new float[x.Length];
            for (var k = 0; k < x.Length; k++)
                sx[k] = (float)((x[k] - mean) / Math.Sqrt(variance));
            return sx;
        }

        private static float[] Divide(double[] x, double divider)
        {
            var sx = new float[x.Length];
            for (var k = 0; k < x.Length; k++)
                sx[k] = (float)(x[k]/divider);
            return sx;
        }


    }
}
