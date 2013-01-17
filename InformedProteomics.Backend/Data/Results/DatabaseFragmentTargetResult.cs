using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.Data.Results
{
	/// <summary>
	/// Contains information related to the resulting search for a specific fragment target.
	/// </summary>
	public class DatabaseFragmentTargetResult
	{
		/// <summary>
		/// The Fragment target used.
		/// </summary>
		public DatabaseFragmentTarget DatabaseFragmentTarget { get; private set; }

		// XY data points representing the elution profile for this target.
		public XYData XYData { get; private set; }

		/// <summary>
		/// Information related to the elution profile for this target.
		/// </summary>
		public XICProfile XICProfile { get; private set; }

		/// <summary>
		/// Contains scores and information as compiled by DeconTools.
		/// </summary>
		public ChromPeakQualityData PeakQualityData { get; private set; }

		public DatabaseFragmentTargetResult(DatabaseFragmentTarget target, XYData xyData, XICProfile xicProfile, ChromPeakQualityData peakQualityData)
		{
			this.DatabaseFragmentTarget = target;
			this.XYData = xyData;
			this.XICProfile = xicProfile;
			this.PeakQualityData = peakQualityData;
		}
	}
}
