using MathNet.Numerics.LinearAlgebra.Single;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Scoring
{
    class FragmentXICSetScorer
    {
        public FragmentXICSet FragmentXICSet { get; private set; }
        public FragmentParameter Parameter { get; private set; }
        public double[] PrecursorXIC { get; private set; }
        public List<IonType> UsedIonTypes { get; private set; }
        public float RawScore { get; private set; }
        public float Score { get; private set; }

        public FragmentXICSetScorer(FragmentXICSet fragmentXICSet, List<IonType> usedIonTypes, double[] precursorXIC, FragmentParameter par)
        {
            FragmentXICSet = fragmentXICSet;
            Parameter = par;
            UsedIonTypes = usedIonTypes;
            PrecursorXIC = precursorXIC;
            Score = GetScore();
        }

        private float GetScore()
        {
            if (FragmentXICSet.Count == 0) return 0;
            var r = GetCorrelationMatrices(UsedIonTypes, Parameter);
            RawScore = new MultipleCorrelationCoefficient(GetX(FragmentXICSet, UsedIonTypes, PrecursorXIC), GetY(PrecursorXIC), r[0], r[1]).Get();
            return SubScoreFactory.GetProductIonXICLikelihoodRatioScore(RawScore, UsedIonTypes.Count, Parameter);
        }

        private static DenseMatrix GetX(FragmentXICSet set, List<IonType> ions, double[] precursorXIC)
        {
            var x = new DenseMatrix(precursorXIC.Length, ions.Count);
            for (var i = 0; i < ions.Count; i++)
            {
                var vs = set.GetIonXIC(ions[i]);
                for (var j=0;j<vs.Length;j++)
                {
                    x.At(j, i, (float)vs[j]);
                }
            }
            return x;
        }

        private static DenseMatrix GetY(double[] precursorXIC)
        {
            var y = new DenseMatrix(precursorXIC.Length, 1);
            for (var i = 0; i < precursorXIC.Length; i++)
            {
                y.At(i, 0, (float)precursorXIC[i]);
            }
                return y;
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
                    rxx.At(i, k, SubScoreFactory.GetProductIonCorrelationCoefficient(ions[i], ions[k], par));
                }
                rxy.At(i, 0, SubScoreFactory.GetProductIonCorrelationCoefficient(ions[i], par));
            }

            return new[] { rxx, rxy };
        }

    }
}
