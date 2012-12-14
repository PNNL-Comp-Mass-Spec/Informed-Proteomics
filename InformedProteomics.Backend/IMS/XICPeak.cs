using DeconTools.Backend.Core;

namespace InformedProteomics.Backend.IMS
{
    public class XICPeak
    {
		/// <summary>
		/// Constructor of XICPeak that uses DeconTools objects to fill in corresponding data.
		/// </summary>
		/// <param name="chromPeak">The DeconTools peak object.</param>
		public XICPeak(ChromPeak chromPeak)
		{
			this.Intensity = chromPeak.Height;
			this.NormalizedElutionTime = chromPeak.NETValue;
		}

        public FrameSet FrameSet
        {
            get
            {
                throw new System.NotImplementedException();
            }
            set
            {
            }
        }

		/// <summary>
		/// Height of the the Peak.
		/// </summary>
    	public float Intensity { get; private set; }

		/// <summary>
		/// Normalized elution time of the Peak.
		/// </summary>
		public double NormalizedElutionTime { get; private set; }

        public XICPeak GetNextXICPeak()
        {
            throw new System.NotImplementedException();
        }

        public XICPeak GetPrevXICPeak()
        {
            throw new System.NotImplementedException();
        }
    }
}
