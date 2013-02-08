using InformedProteomics.Backend.Data.Results;
using MathNet.Numerics.LinearAlgebra.Single;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Scoring
{
    class PrecursorIonScorer
    {
        public List<DatabaseSubTargetResult> PrecursorResults { get; private set; }
        public DatabaseSubTargetResult PrecursorResultRep { get; private set; }
        public List<int> ChargeStateList { get; private set; }
        public float Score { get; private set; }

        public PrecursorIonScorer(DatabaseMultipleSubTargetResult matchedResult)
        {
            PrecursorResults = matchedResult.SubTargetResultList;
            PrecursorResultRep = matchedResult.PrecursorResultRep;
            ChargeStateList = matchedResult.ChargeStateList;
            Score = GetScore();
            //Console.WriteLine("Precursor Score : " + Score);
        }


        private float GetScore()
        {
            if (PrecursorResults.Count == 1) return 0;

            var r = GetCorrelationMatrices();
            var rawScore = new MultipleCorrelationCoefficient(GetX(), GetY(), r[0], r[1]).Get();
            return SubScoreFactory.GetPrecursorIonLikelihoodRatioScore(rawScore, PrecursorResults.Count);
        }

        private DenseMatrix GetX()
        {
            var j = PrecursorResults.Count;
            var l = PrecursorResultRep.XYData.Yvalues.Length;

            var x = new DenseMatrix(l, j-1, 0);

            var n = 0;
            for (var i = 0; i < j;i++)
            {
                var r = PrecursorResults.ElementAt(i);
                if (r.Equals(PrecursorResultRep)) continue;

                var s = r.XYData.Yvalues;

                for (var k = 0; k < l; k++)
                {
                    x.At(k,n,(float)s[k]);
                }
                n++;
            }

            return x;
        }

        private DenseMatrix GetY()
        {
            var l = PrecursorResultRep.XYData.Yvalues.Length;

            var y = new DenseMatrix(l, 1, 0);
            var s = PrecursorResultRep.XYData.Yvalues;

            for (var j = 0; j < l; j++)
            {
                y.At(j, 0, (float)s[j]);
            }
             
            return y;
        }

        private DenseMatrix[] GetCorrelationMatrices()
        {
            var j = PrecursorResults.Count;
            var rxx = new DenseMatrix(j - 1, j - 1, 0);
            var rxy = new DenseMatrix(j - 1, 1, 0);

            var charges = new int[j-1];
            var charge = 0;
            var n = 0;
            for (var i = 0; i < j; i++)
            {
                var r = PrecursorResults.ElementAt(i);
                if (r.Equals(PrecursorResultRep))
                {
                    charge = ChargeStateList.ElementAt(i);
                    continue;
                }
                charges[n] = ChargeStateList.ElementAt(i);
                n++;
            }

            for (var i = 0; i < j - 1;i++)
            {
                for (var k = 0; k < j - 1; k++)
                {
                    rxx.At(i,k,SubScoreFactory.GetPrecursorIonCorrelationCoefficient(charges[i], charges[k]));
                }
                rxy.At(i, 0, SubScoreFactory.GetPrecursorIonCorrelationCoefficient(charges[i], charge));
            }

            return new[] { rxx, rxy }; 
        }

      
        

       
    }
}
