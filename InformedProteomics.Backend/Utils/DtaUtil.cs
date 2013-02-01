using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.Utils
{
	public class DtaUtil
	{
		public static void AppendToDtaFile(DatabaseMultipleSubTargetResult result, TextWriter dtaWriter)
		{
			string dtaString = CreateSingleDtaEntry(result);

			// Only write out valid strings
			if(dtaString != null && !dtaString.Equals(""))
			{
				dtaWriter.WriteLine(dtaString);
			}

			// Flush the DTA so it can be read immediately
			dtaWriter.Flush();
		}

		private static string CreateSingleDtaEntry(DatabaseMultipleSubTargetResult result)
		{
			DatabaseSubTargetResult precursorRep = result.PrecursorResultRep;
			DatabaseSubTarget precursorTarget = precursorRep.DatabaseSubTarget;
			XICProfile xicProfile = precursorRep.XICProfile;

			string sequence = precursorTarget.Code;
			string empiricalFormula = precursorTarget.EmpiricalFormula;
			double elutionTime = Math.Round(xicProfile.ApexPeak.NormalizedElutionTime, 4);
			double abundance = xicProfile.Abundance;
			int chargeState = precursorTarget.ChargeState;
			double mass = precursorTarget.MonoIsotopicMass;

			StringBuilder dtaString = new StringBuilder();
			String header = "=================================== \"" + sequence + "." + empiricalFormula + "." + chargeState + "." + elutionTime + "." + abundance + ".dta\" ==================================";
			String parentIonLine = (mass + Globals.Hydrogen_MASS) + " " + chargeState;

			dtaString.AppendLine(header);
			dtaString.AppendLine(parentIonLine);

			List<MSPeak> msPeakList = new List<MSPeak>();

			foreach (DatabaseFragmentTargetResult fragmentResult in result.FragmentResultList)
			{
				if (fragmentResult == null || fragmentResult.PeakQualityData == null || fragmentResult.PeakQualityData.IsotopicProfile == null || fragmentResult.PeakQualityData.IsotopicProfile.Peaklist == null) continue;
				msPeakList.AddRange(fragmentResult.PeakQualityData.IsotopicProfile.Peaklist);
			}

			// If no peaks, don't return the dta string
			if(!msPeakList.Any())
			{
				return "";
			}

			double previousXValue = 0;
			foreach (var msPeak in msPeakList.OrderBy(x => x.XValue))
			{
				if (Math.Abs(msPeak.XValue - previousXValue) < 0.001) continue;
				previousXValue = msPeak.XValue;

				dtaString.AppendLine(msPeak.XValue + " " + msPeak.Height);
			}

			return dtaString.ToString();
		}
	}
}
