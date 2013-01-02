using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NETPrediction;

namespace InformedProteomics.Backend.Utils
{
	public class PeptideUtil
	{
		private static readonly ElutionTimePredictionKangas _netPrediction;

		static PeptideUtil()
		{
			_netPrediction = new ElutionTimePredictionKangas();
		}

		/// <summary>
		/// Uses the pnet calculator to get a predicted normalized elution time for a given peptide sequence.
		/// </summary>
		/// <param name="peptideString">The peptide string can either be in the form PEPTIDE or X.PEPTIDE.X</param>
		/// <returns>The predicted normalized elution time value.</returns>
		public static float CalculatePredictedElutionTime(String peptideString)
		{
			return _netPrediction.GetElutionTime(PrepareSequenceForNetPrediction(peptideString));
		}

		private static string PrepareSequenceForNetPrediction(string peptideSequence)
		{
			return peptideSequence.Replace('C', 'c');
		}
	}
}
