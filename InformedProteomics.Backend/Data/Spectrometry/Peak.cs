using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Peak: IComparable<Peak>
    {
        public Peak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public double Mz { get; private set; }
        public double Intensity { get; private set; }

        public int CompareTo(Peak other)
        {
	        return Mz.CompareTo(other.Mz);
        }

    	public bool Equals(Peak other)
    	{
    		if (ReferenceEquals(null, other)) return false;
    		if (ReferenceEquals(this, other)) return true;
			return Math.Abs(Mz - other.Mz) < 1e-9;
    	}

    	public override bool Equals(object obj)
    	{
    		if (ReferenceEquals(null, obj)) return false;
    		if (ReferenceEquals(this, obj)) return true;
    		if (obj.GetType() != typeof (Peak)) return false;
    		return Equals((Peak) obj);
    	}

    	public override int GetHashCode()
    	{
    		return Mz.GetHashCode();
    	}

    	public static bool operator ==(Peak left, Peak right)
    	{
    		return Equals(left, right);
    	}

    	public static bool operator !=(Peak left, Peak right)
    	{
    		return !Equals(left, right);
    	}
    }

    /// <summary>
    /// Sort by reverse order of intensities (highest intensity comes first)
    /// </summary>
    public class IntensityComparer : IComparer<Peak>
    {
        public int Compare(Peak x, Peak y)
        {
	        return y.Intensity.CompareTo(x.Intensity);
        }
    }

    /// <summary>
    /// Compare by m/z. Two peaks within ppmTolerance are considered to be equal.
    /// </summary>
    public class MzComparerWithPpmTolerance : IComparer<Peak>, IEqualityComparer<Peak>
    {
        public MzComparerWithPpmTolerance(double ppmTolerance)
        {
            _ppmTolerance = ppmTolerance;
        }

        public int Compare(Peak x, Peak y)
        {
            if (Equals(x, y)) return 0;
            return x.Mz.CompareTo(y.Mz);
        }

        public bool Equals(Peak x, Peak y)
        {
            return Math.Abs((x.Mz - y.Mz)/x.Mz*1E6) <= _ppmTolerance;
        }

        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly double _ppmTolerance;
    }

}
