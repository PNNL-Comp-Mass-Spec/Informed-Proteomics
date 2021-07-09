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
        public int ScanNum { get; }

        /// <summary>
        /// Compare two LcMsPeaks (for sorting)
        /// </summary>
        /// <param name="other"></param>
        /// <returns>0, -1 or 1</returns>
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
