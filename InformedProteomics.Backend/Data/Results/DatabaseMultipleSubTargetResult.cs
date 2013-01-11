using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Results
{
	public class DatabaseMultipleSubTargetResult
	{
		public double ElutionTime { get; private set; }
		public List<DatabaseSubTargetResult> SubTargetResultList { get; private set; }
		public List<int> ChargeStateList { get; private set; }
		public DatabaseSubTargetResult PrecursorResultRep { get; private set; }
		public IList<DatabaseFragmentTargetResult> FragmentResultList { get; set; } 

		public DatabaseMultipleSubTargetResult(DatabaseSubTargetResult result)
		{
			this.ElutionTime = result.XICProfile.ApexPeak.NormalizedElutionTime;
			this.SubTargetResultList = new List<DatabaseSubTargetResult> { result };
			this.ChargeStateList = new List<int> { result.DatabaseSubTarget.ChargeState };
			this.PrecursorResultRep = result;
			this.FragmentResultList = new List<DatabaseFragmentTargetResult>();
		}

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

		public void AddNewResult(DatabaseSubTargetResult result)
		{
			this.SubTargetResultList.Add(result);
			this.ChargeStateList.Add(result.DatabaseSubTarget.ChargeState);
			this.ElutionTime = CalculateAverageElutionTime();

			if(result.XICProfile.ApexPeak.Intensity > this.PrecursorResultRep.XICProfile.ApexPeak.Intensity)
			{
				this.PrecursorResultRep = result;
			}
		}

		private double CalculateAverageElutionTime()
		{
			return SubTargetResultList.Select(databaseSubTargetResult => databaseSubTargetResult.XICProfile.ApexPeak.NormalizedElutionTime).Average();
		}
	}
}
