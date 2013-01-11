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
	public class DatabaseFragmentTargetResult
	{
		public DatabaseFragmentTarget DatabaseFragmentTarget { get; private set; }
		public XYData XYData { get; private set; }
		public XICProfile XICProfile { get; private set; }
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
