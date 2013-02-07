using System.Linq;
using MathNet.Numerics.LinearAlgebra.Generic;
using MathNet.Numerics.LinearAlgebra.Single;
using System;

namespace InformedProteomics.Backend.Scoring
{
    public class MultipleCorrelationCoefficient
    {
        private DenseMatrix X { get; set; }
        private DenseMatrix Y { get; set; }
        private DenseMatrix Rxx { get; set; }
        private DenseMatrix Rxy { get; set; }

        public MultipleCorrelationCoefficient(DenseMatrix x, DenseMatrix y, DenseMatrix rxx, DenseMatrix rxy)// x should be appended to 1s
        {
            X = Standardize(x);
            Y = Standardize(y);
            Rxx = rxx;
            Rxy = rxy;
        }

        private void PrintMatrix(DenseMatrix m)
        {
            for (var i = 0; i < m.RowCount; i++)
            {
                for (var j = 0; j < m.ColumnCount; j++)
                {
                    Console.Write(m.At(i,j) + "\t");
                }
                Console.WriteLine();
            }
        }

        public float Get()
        {
            var n = Y.RowCount;
           // var j = X.ColumnCount;
           
            var estY = X.Multiply(Rxx.Inverse().Multiply(Rxy));
            var yvar = new DenseMatrix(n, 1, 1).Transpose().Multiply(Y).At(0,0);
            yvar = yvar * yvar / n;

            var ssr = estY.Transpose().Multiply(Y).At(0, 0) - yvar;
            var sst = Y.Transpose().Multiply(Y).At(0, 0) - yvar;
           
            
            var r2 = Math.Max(0,ssr) / sst;
            
            return r2;
        }

        private static DenseMatrix Standardize(Matrix<float> x) // return standardized x (i.e., mean = 0, var = 1) 
        {
            var sx = new DenseMatrix(x.RowCount, x.ColumnCount);
            for (var i = 0; i < x.ColumnCount; i++)
            { 
                var c = x.Column(i);
                var m = GetSampleMean(c);
                var v = GetSampleVariance(c, m);
               
                for (var k = 0; k < x.RowCount; k++)
                {
                    sx.At(k,i, (x.At(k,i) - m) / (float)Math.Sqrt(v));
                }
            }
            return sx;
        }


        private static float GetSampleMean(Vector<float> x)
        {
            var m = Enumerable.Sum(x);
            return m / x.Count;
        }

        private static float GetSampleVariance(Vector<float> x, float m)
        {
            var var = x.Sum(v => (v - m)*(v - m));
            return var / (x.Count - 1);
        }

        /*
        static public void Main(String[] args)
        {
            var x = new DenseMatrix(18, 2, new float[] { 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 2, 2, 4, 4, 8, 8, 2, 2, 4, 4, 8, 8, 2, 2, 4, 4, 8, 8});//{ 2, 2, 2, 2, 2, 4, 2, 4, 2, 8, 2, 8, 4, 2, 4, 2, 4, 4, 4, 4, 4, 8, 4, 8, 8, 2, 8, 2, 8, 4, 8, 4, 8, 8, 8, 8 });
            var y = new DenseMatrix(18, 1, new float[] { 35, 39, 21, 31, 6, 8, 40, 52, 34, 42, 18, 26, 61, 73, 58, 66, 46, 52 });
            var n = y.RowCount;
            var j = x.ColumnCount;
            var ones = new DenseMatrix(n, 1, 1);
            var newX = new DenseMatrix(n, j + 1, 0);
            ones.Append(x, newX);

            var rxx = (DenseMatrix)newX.TransposeThisAndMultiply(newX);
            var rxy = (DenseMatrix)newX.TransposeThisAndMultiply(y);
            var r2 = new MultipleCorrelationCoefficient(newX, y, rxx, rxy).Get();
            Console.Write(r2);
        }
        */
    }
}
