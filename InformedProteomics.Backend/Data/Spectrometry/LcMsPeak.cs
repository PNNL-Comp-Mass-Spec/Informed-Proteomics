using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// LcMsPeak: like <see cref="Peak"/>, but with a scan number
    /// </summary>
    public class LcMsPeak : Peak, IComparable<LcMsPeak>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="intensity"></param>
        /// <param name="scanNum"></param>
        public LcMsPeak(double mz, double intensity, int scanNum)
            : base(mz, intensity)
        {
            ScanNum = scanNum;
        }

        /// <summary>
        /// Scan number where peak was observed
        /// </summary>
        public int ScanNum { get; private set; }

        /// <summary>
        /// Replace the data in the peak (used by PbfLcMsRun)
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="intensity"></param>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        public LcMsPeak ReplaceData(double mz, double intensity, int scanNum)
        {
            Mz = mz;
            Intensity = intensity;
            ScanNum = scanNum;
            return this;
        }

        /// <summary>
        /// Compare 2 LcMsPeaks (for sorting)
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(LcMsPeak other)
        {
            var mzCompare = Mz.CompareTo(other.Mz);
            if (mzCompare == 0)
            {
                // Force a stable sort
                return ScanNum.CompareTo(other.ScanNum);
            }
            return mzCompare;
        }
    }
}
