using System.Collections.Generic;
using DeconTools.Backend.Core;

namespace InformedProteomics.Backend.IMS
{
    public class XICProfile
    {
		/// <summary>
		/// A list of all peaks contained in the profile.
		/// </summary>
		public IList<XICPeak> XICPeakList { get; private set; }

		/// <summary>
		/// Link to the peak data for the apex of the profile.
		/// </summary>
		public XICPeak ApexPeak { get; private set; }

		/// <summary>
		/// Isotopic Fit Score as defined by DeconTools. Dot product of observed isotopic profile compared to theoretical isotopic profile.
		/// </summary>
		public double DeconToolsFitScore { get; private set; }

		/// <summary>
		/// Area under the curve.
		/// </summary>
		public double Abundance { get; private set; }

		/// <summary>
		/// Score defined by DeconTools that measures the level of intensity of nearby peaks that could possibly interfere with the peaks of the profile.
		/// </summary>
		public double InterferenceScore { get; private set; }

		/// <summary>
		/// Constructor of XICProfile that uses DeconTools objects to fill in corresponding data.
		/// </summary>
		/// <param name="peakQualityData">This object stores various scores of the XICProfile as defined by DeconTools.</param>
		/// <param name="apexPeak">The apex of the profile as found by DeconTools.</param>
		public XICProfile(ChromPeakQualityData peakQualityData, ChromPeak apexPeak)
		{
			this.DeconToolsFitScore = peakQualityData.FitScore;
			this.Abundance = peakQualityData.Abundance;
			this.InterferenceScore = peakQualityData.InterferenceScore;

			this.ApexPeak = new XICPeak(apexPeak);
			this.XICPeakList = new List<XICPeak> { this.ApexPeak };
		}

        public XICPeak GetApexXICPeak()
        {
        	return this.ApexPeak;
        }
    }
}
