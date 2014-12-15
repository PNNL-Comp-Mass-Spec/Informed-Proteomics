using System;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;

namespace InformedProteomics.Backend.Utils
{
	public class SavitzkyGolaySmoother
	{
		private readonly int _mNumPointsForSmoothing;
		private readonly Matrix<double> _mSmoothingFiltersConjugateTranspose;

		public SavitzkyGolaySmoother(int pointsForSmoothing, int polynomialOrder)
		{
			if (pointsForSmoothing < 3) throw new ArgumentOutOfRangeException("savGolayPoints must be an odd number 3 or higher");
			if (pointsForSmoothing % 2 == 0) throw new ArgumentOutOfRangeException("savGolayPoints must be an odd number 3 or higher");

			_mNumPointsForSmoothing = pointsForSmoothing;
			var smoothingFilters = CalculateSmoothingFilters(polynomialOrder, pointsForSmoothing);
			_mSmoothingFiltersConjugateTranspose = smoothingFilters.ConjugateTranspose();
		}

		public void Smooth(ref double[,] inputValues)
		{
			// TODO: Using the matrix works, but does a lot of data accesses. Can improve by working out all the data access myself? I might be able to cut down on number of data accesses, but not sure.
		    var inputMatrix = new DenseMatrix(inputValues.GetLength(0), inputValues.GetLength(1));

			for (int i = 0; i < inputMatrix.RowCount; i++)
			{
				inputMatrix.SetRow(i, Smooth(inputMatrix.Row(i).ToArray()));
			}

			for (int i = 0; i < inputMatrix.ColumnCount; i++)
			{
				inputMatrix.SetColumn(i, Smooth(inputMatrix.Column(i).ToArray()));
			}

			inputValues = inputMatrix.ToArray();
		}

		public double[] Smooth(double[] inputValues)
		{
			// No need to smooth if all values are 0
			if (IsEmpty(inputValues)) return inputValues;

			int m = (_mNumPointsForSmoothing - 1) / 2;
			int colCount = inputValues.Length;
			double[] returnYValues = new double[colCount];

			var conjTransposeMatrix = _mSmoothingFiltersConjugateTranspose;

			for (int i = 0; i <= m; i++)
			{
				var conjTransposeColumn = conjTransposeMatrix.Column(i);

				double multiplicationResult = 0;
				for (int z = 0; z < _mNumPointsForSmoothing; z++)
				{
					multiplicationResult += (conjTransposeColumn[z] * inputValues[z]);
				}

				returnYValues[i] = multiplicationResult;
			}

			var conjTransposeColumnResult = conjTransposeMatrix.Column(m);

			for (int i = m + 1; i < colCount - m - 1; i++)
			{
				double multiplicationResult = 0;
				for (int z = 0; z < _mNumPointsForSmoothing; z++)
				{
					multiplicationResult += (conjTransposeColumnResult[z] * inputValues[i - m + z]);
				}
				returnYValues[i] = multiplicationResult;
			}

			for (int i = 0; i <= m; i++)
			{
				var conjTransposeColumn = conjTransposeMatrix.Column(m + i);

				double multiplicationResult = 0;
				for (int z = 0; z < _mNumPointsForSmoothing; z++)
				{
					multiplicationResult += (conjTransposeColumn[z] * inputValues[colCount - _mNumPointsForSmoothing + z]);
				}
				returnYValues[colCount - m - 1 + i] = multiplicationResult;
			}

			return returnYValues;
		}

		private DenseMatrix CalculateSmoothingFilters(int polynomialOrder, int filterLength)
		{
			int m = (filterLength - 1) / 2;
			var denseMatrix = new DenseMatrix(filterLength, polynomialOrder + 1);

			for (int i = -m; i <= m; i++)
			{
				for (int j = 0; j <= polynomialOrder; j++)
				{
					denseMatrix[i + m, j] = Math.Pow(i, j);
				}
			}

			var sTranspose = (DenseMatrix)denseMatrix.ConjugateTranspose();
			var f = sTranspose * denseMatrix;
			var fInverse = (DenseMatrix)f.LU().Solve(DenseMatrix.CreateIdentity(f.ColumnCount));
			var smoothingFilters = denseMatrix * fInverse * sTranspose;

			return smoothingFilters;
		}

		private bool IsEmpty(double[] inputValues)
		{
			foreach (var inputValue in inputValues)
			{
				if (inputValue > 0) return false;
			}

			return true;
		}
	}
}
