using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.Data.Results
{
	/// <summary>
	/// Contains information related to the resulting search for a single target.
	/// </summary>
	public class DatabaseSubTargetResult
	{
		/// <summary>
		/// The original target, not including the specific charge state.
		/// </summary>
		public DatabaseTarget DatabaseTarget { get; private set; }

		/// <summary>
		/// The actual target that was used.
		/// </summary>
		public DatabaseSubTarget DatabaseSubTarget { get; private set; }

		/// <summary>
		/// The isotopic fit score, as calculated by DeconTools
		/// </summary>
		public double IsotopicFitScore { get; private set; }

		/// <summary>
		/// XY data points that repesent the elution profile for this target.
		/// </summary>
		public XYData XYData { get; private set; }

		/// <summary>
		/// Information related to the elution profile for this target.
		/// </summary>
		public XICProfile XICProfile { get; private set; }

		public DatabaseSubTargetResult(DatabaseSubTarget subTarget, DatabaseTarget databaseTarget, XYData xyData, XICProfile xicProfile, double isotopicFitScore)
		{
			this.DatabaseSubTarget = subTarget;
			this.DatabaseTarget = databaseTarget;
			this.XYData = xyData;
			this.XICProfile = xicProfile;
			this.IsotopicFitScore = isotopicFitScore;
		}
	}
}
