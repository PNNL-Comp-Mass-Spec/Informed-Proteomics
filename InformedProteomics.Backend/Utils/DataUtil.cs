using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Algorithms;

namespace InformedProteomics.Backend.Utils
{
	public class DataUtil
	{
		private static ChromatogramCorrelator chromatogramCorrelator;

		static DataUtil()
		{
			chromatogramCorrelator = new ChromatogramCorrelator();
		}

		public static void CorrelateXYData(XYData profile1, XYData profile2, int diffBetweenXValues, out double slope, out double intercept, out double rSquared)
		{
			XYData alignedProfile = AlignXYData(profile1, profile2, diffBetweenXValues);
			chromatogramCorrelator.GetElutionCorrelationData(profile1, alignedProfile, out slope, out intercept, out rSquared);
		}

		public static XYData AlignXYData(XYData referenceProfile, XYData profileToAlign, int diffBetweenXValues)
		{
			double[] referenceXValues = referenceProfile.Xvalues;
			double[] oldYValues = profileToAlign.Yvalues;
			double[] newYValues = new double[referenceXValues.Length];
			List<double> oldXValues = profileToAlign.Xvalues.ToList();

			for(int i = 0; i < referenceXValues.Length; i++)
			{
				int binarySearchResult = oldXValues.BinarySearch(referenceXValues[i] + diffBetweenXValues);

				if(binarySearchResult >= 0)
				{
					newYValues[i] = oldYValues[binarySearchResult];
				}
				else
				{
					newYValues[i] = 0;
				}
			}

			XYData newXYData = new XYData();
			newXYData.SetXYValues(ref referenceXValues, ref newYValues);

			return newXYData;
		}
	}
}
