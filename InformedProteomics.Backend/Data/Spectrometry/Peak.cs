using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using MathNet.Numerics.Distributions;

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
    public class MzComparerWithTolerance : IComparer<Peak>, IEqualityComparer<Peak>
    {
        public MzComparerWithTolerance(Tolerance tolerance)
        {
            if (tolerance.GetUnit() == ToleranceUnit.Ppm)
            {
                _equalityComparer = new MzComparerWithPpmTolerance(tolerance.GetValue());
            }
            else if (tolerance.GetUnit() == ToleranceUnit.Th)
            {
                _equalityComparer = new MzComparerWithToleranceMz(tolerance.GetValue());
            }
            else throw new NotSupportedException("The tolerance unite must be ppm or Th.");
        }

        public MzComparerWithTolerance(double ppmTolerance): this(new Tolerance(ppmTolerance))
        {
        }

        public int Compare(Peak x, Peak y)
        {
            if (Equals(x, y)) return 0;
            return x.Mz.CompareTo(y.Mz);
        }

        public bool Equals(Peak x, Peak y)
        {
            return _equalityComparer.Equals(x, y);
        }

        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly IEqualityComparer<Peak> _equalityComparer;
    }

    /// <summary>
    /// Compare by m/z. Two peaks within ppmTolerance are considered to be equal.
    /// </summary>
    public class MzComparerWithPpmTolerance : IEqualityComparer<Peak>
    {
        public MzComparerWithPpmTolerance(double ppmTolerance)
        {
            _ppmTolerance = ppmTolerance;
        }

        public bool Equals(Peak x, Peak y)
        {
            return Math.Abs((x.Mz - y.Mz) / x.Mz * 1E6) <= _ppmTolerance;
        }

        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly double _ppmTolerance;
    }

    /// <summary>
    /// Compare by m/z. Two peaks within toleranceTh Th are considered to be equal.
    /// </summary>
    public class MzComparerWithToleranceMz : IEqualityComparer<Peak>
    {
        public MzComparerWithToleranceMz(double toleranceTh)
        {
            _toleranceTh = toleranceTh;
        }

        public bool Equals(Peak x, Peak y)
        {
            return Math.Abs((x.Mz - y.Mz)) <= _toleranceTh;
        }

        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly double _toleranceTh;
    }

    /// <summary>
    /// Compare by m/z. Two peaks within ppmTolerance are considered to be equal.
    /// </summary>
    public class MzComparerWithBinning : IComparer<Peak>, IEqualityComparer<Peak>
    {
        // 27 bits: max error = 16 ppm, 28 bits (8 ppm), 26 bits (32 ppm)
        public MzComparerWithBinning(int numBits = 27)
        {
            _numShifts = sizeof(double)*8 - numBits;
        }

        public bool Equals(Peak x, Peak y)
        {
            return GetBinNumber(x.Mz) == GetBinNumber(y.Mz);
        }

        public bool Equals(double x, double y)
        {
            return GetBinNumber(x) == GetBinNumber(y);
        }

        public int GetHashCode(Peak p)
        {
            return GetBinNumber(p.Mz);
        }

        public int Compare(Peak x, Peak y)
        {
            return GetBinNumber(x.Mz).CompareTo(GetBinNumber(y.Mz));
        }

        public double GetRoundedValue(double mz)
        {
            var converted = BitConverter.DoubleToInt64Bits(mz);
            var rounded = (converted >> _numShifts) << _numShifts;
            var roundedDouble = BitConverter.Int64BitsToDouble(rounded);
            return roundedDouble;
        }

        public int GetBinNumber(double mz)
        {
            var converted = BitConverter.DoubleToInt64Bits(mz);
            return (int) (converted >> _numShifts);
            //var rounded = (converted >> _numShifts) << _numShifts;
            //return (int)(rounded >> _numShifts);
        }

        // inclusive
        public double GetMzStart(int binNum)
        {
            var rounded = (long)binNum << _numShifts;
            return BitConverter.Int64BitsToDouble(rounded);
        }

        // exclusive
        public double GetMzEnd(int binNum)
        {
            var rounded = (long)(binNum+1) << _numShifts;
            return BitConverter.Int64BitsToDouble(rounded);
        }

        public double GetMzAverage(int binNum)
        {
            return (GetMzStart(binNum) + GetMzEnd(binNum))*0.5;
        }

        public Tolerance GetTolerance()
        {
            // 27 bits: max error = 16 ppm, 28 bits (8 ppm), 26 bits (32 ppm)
            var numBits = sizeof(double) * 8 - _numShifts;
            var ppm = 16*Math.Pow(2, 27 - numBits);

            return new Tolerance(ppm);
        }

        private readonly int _numShifts;
    }

}
