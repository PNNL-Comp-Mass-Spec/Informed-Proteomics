using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.Data.Results
{
	public class DatabaseSubTargetResult
	{
		public DatabaseTarget DatabaseTarget { get; private set; }
		public DatabaseSubTarget DatabaseSubTarget { get; private set; }
		public XYData XYData { get; private set; }
		public XICProfile XICProfile { get; private set; }
		public IList<DatabaseFragmentTargetResult> FragmentResultList { get; set; } 

		public DatabaseSubTargetResult(DatabaseSubTarget subTarget, DatabaseTarget databaseTarget, XYData xyData, XICProfile xicProfile)
		{
			this.DatabaseSubTarget = subTarget;
			this.DatabaseTarget = databaseTarget;
			this.XYData = xyData;
			this.XICProfile = xicProfile;
			this.FragmentResultList = new List<DatabaseFragmentTargetResult>();
		}
	}
}
