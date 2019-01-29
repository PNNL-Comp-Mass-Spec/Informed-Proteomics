using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// A peak in an XIC
    /// </summary>
    public class XicPoint: IComparable<XicPoint>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="mz"></param>
        /// <param name="intensity"></param>
        public XicPoint(int scanNum, double mz, double intensity)
        {
            ScanNum = scanNum;
            Mz = mz;
            Intensity = intensity;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peak"></param>
        public XicPoint(LcMsPeak peak)
        {
            ScanNum = peak.ScanNum;
            Mz = peak.Mz;
            Intensity = peak.Intensity;
        }

        /// <summary>
        /// Peak scan number
        /// </summary>
        public int ScanNum { get; }

        /// <summary>
        /// Peak m/z
        /// </summary>
        public double Mz { get; }

        /// <summary>
        /// Peak intensity
        /// </summary>
        public double Intensity { get; }

        /// <summary>
        /// Compare 2 XicPoints
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(XicPoint other)
        {
            return ScanNum - other.ScanNum;
        }

        /// <summary>
        /// Check 2 XicPoints for equality
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(XicPoint other)
        {
            if (other == null) return false;
            if (ReferenceEquals(this, other)) return true;
            return other.ScanNum == ScanNum && Math.Abs(other.Mz - Mz) < 0.001 && Math.Abs(other.Intensity - Intensity) < 0.001;
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return ScanNum + "," + Mz + "," + Intensity;
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (obj == null) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != typeof(XicPoint)) return false;
            return Equals((XicPoint)obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return ScanNum;
        }

        /// <summary>
        /// Overloaded equality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator ==(XicPoint left, XicPoint right)
        {
            return Equals(left, right);
        }

        /// <summary>
        /// Overloaded inequality operator
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator !=(XicPoint left, XicPoint right)
        {
            return !Equals(left, right);
        }
    }
}
