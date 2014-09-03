using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class XicPoint: IComparable<XicPoint>
    {
        public XicPoint(int scanNum, double mz, double intensity)
        {
            ScanNum = scanNum;
            Mz = mz;
            Intensity = intensity;
        }

        public int ScanNum { get; private set; }
        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public int CompareTo(XicPoint other)
        {
            return ScanNum - other.ScanNum;
        }

		public bool Equals(XicPoint other)
		{
			if (ReferenceEquals(null, other)) return false;
			if (ReferenceEquals(this, other)) return true;
			return other.ScanNum == ScanNum && Math.Abs(other.Mz - Mz) < 0.001 && Math.Abs(other.Intensity - Intensity) < 0.001;
		}

        public override string ToString()
        {
            return ScanNum + "," + Mz + "," + Intensity;
        }

        public override bool Equals(object obj)
		{
			if (ReferenceEquals(null, obj)) return false;
			if (ReferenceEquals(this, obj)) return true;
			if (obj.GetType() != typeof(XicPoint)) return false;
			return Equals((XicPoint)obj);
		}

		public override int GetHashCode()
		{
			return ScanNum;
		}

		public static bool operator ==(XicPoint left, XicPoint right)
		{
			return Equals(left, right);
		}

		public static bool operator !=(XicPoint left, XicPoint right)
		{
			return !Equals(left, right);
		}
    }
}
