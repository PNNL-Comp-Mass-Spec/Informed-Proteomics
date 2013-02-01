using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Results
{
	/// <summary>
	/// Contains all the results for a peptide target at a specific elution time. 
	/// If the peptide target was found in multiple elution times, then there will be more than 1 of these objects created.
	/// </summary>
	public class DatabaseMultipleSubTargetResult
	{
		/// <summary>
		/// The specific elution time of this result
		/// </summary>
		public double ElutionTime { get; private set; }

		/// <summary>
		/// The best isotopic fit score for this result, as calculated by DeconTools
		/// </summary>
		public double IsotopicFitScore { get; private set; }

		/// <summary>
		/// Each Precursor target result that is attached to this result. Should be 1 per charge state found.
		/// </summary>
		public List<DatabaseSubTargetResult> SubTargetResultList { get; private set; }

		/// <summary>
		/// List of charge states that were discovered when looking for the precursor targets.
		/// </summary>
		public List<int> ChargeStateList { get; private set; }

		/// <summary>
		/// The respresentative precursor target result. This is the most intense precursor discovered.
		/// </summary>
		public DatabaseSubTargetResult PrecursorResultRep { get; private set; }

		/// <summary>
		/// List of Fragment results.
		/// </summary>
		public IList<DatabaseFragmentTargetResult> FragmentResultList { get; set; } 

		public DatabaseMultipleSubTargetResult(DatabaseSubTargetResult result)
		{
			this.ElutionTime = result.XICProfile.ApexPeak.NormalizedElutionTime;
			this.SubTargetResultList = new List<DatabaseSubTargetResult> { result };
			this.ChargeStateList = new List<int> { result.DatabaseSubTarget.ChargeState };
			this.PrecursorResultRep = result;
			this.FragmentResultList = new List<DatabaseFragmentTargetResult>();
			this.IsotopicFitScore = result.IsotopicFitScore;
		}

		/// <summary>
		/// Checks to see if a new precursor target result fits in with this result object.
		/// </summary>
		/// <param name="result">The precursor target result to test.</param>
		/// <returns>True if the new result belongs, false otherwise.</returns>
		public bool DoesNewResultBelong(DatabaseSubTargetResult result)
		{
			double elutionTime = result.XICProfile.ApexPeak.NormalizedElutionTime;

			// Make sure elution time is close enough
			if(Math.Abs(this.ElutionTime - elutionTime) < 0.005)
			{
				int chargeState = result.DatabaseSubTarget.ChargeState;

				// Make sure we have not included this charge state yet
				if (!this.ChargeStateList.Contains(chargeState))
				{
					// TODO: Also correlate the peak shape?
					return true;
				}
			}

			return false;
		}

		/// <summary>
		/// Adds a new precursor result to this result. This method will update the list of charge states, average elution time, and the representative precursor.
		/// </summary>
		/// <param name="result"></param>
		public void AddNewResult(DatabaseSubTargetResult result)
		{
			this.SubTargetResultList.Add(result);
			this.ChargeStateList.Add(result.DatabaseSubTarget.ChargeState);
			this.ElutionTime = CalculateAverageElutionTime();

			if(result.XICProfile.ApexPeak.Intensity > this.PrecursorResultRep.XICProfile.ApexPeak.Intensity)
			{
				this.PrecursorResultRep = result;
			}

			if(result.IsotopicFitScore < this.IsotopicFitScore)
			{
				this.IsotopicFitScore = result.IsotopicFitScore;
			}
		}

		private double CalculateAverageElutionTime()
		{
			return SubTargetResultList.Select(databaseSubTargetResult => databaseSubTargetResult.XICProfile.ApexPeak.NormalizedElutionTime).Average();
		}
	}
}
