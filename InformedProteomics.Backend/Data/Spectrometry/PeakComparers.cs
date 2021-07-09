using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Sort by reverse order of intensities (highest intensity comes first)
    /// </summary>
    public class IntensityComparer : IComparer<Peak>
    {
        /// <summary>
        /// Compare two peaks by intensity
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>0, 1, or -1</returns>
        public int Compare(Peak x, Peak y)
        {
            if (x == null)
            {
                return y == null ? 0 : -1;
            }

            if (y == null)
            {
                return 1;
            }

            return y.Intensity.CompareTo(x.Intensity);
        }
    }

    /// <summary>
    /// Compare by m/z; two peaks within ppmTolerance are considered to be equal
    /// </summary>
    public class MzComparerWithTolerance : IComparer<Peak>, IEqualityComparer<Peak>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="tolerance"></param>
        public MzComparerWithTolerance(Tolerance tolerance)
        {
            if (tolerance.GetUnit() == ToleranceUnit.Ppm)
            {
                _equalityComparer = new MzComparerWithPpmTolerance(tolerance.GetValue());
            }
            else if (tolerance.GetUnit() == ToleranceUnit.Mz)
            {
                _equalityComparer = new MzComparerWithToleranceMz(tolerance.GetValue());
            }
            else
            {
                throw new NotSupportedException("The tolerance unite must be ppm or Th.");
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ppmTolerance"></param>
        public MzComparerWithTolerance(double ppmTolerance) : this(new Tolerance(ppmTolerance))
        {
        }

        /// <summary>
        /// Compare two peaks by m/z with a tolerance
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>0, 1, or-1</returns>
        public int Compare(Peak x, Peak y)
        {
            if (x == null)
            {
                return y == null ? 0 : -1;
            }

            if (y == null)
            {
                return 1;
            }

            if (Equals(x, y))
            {
                return 0;
            }

            return x.Mz.CompareTo(y.Mz);
        }

        /// <summary>
        /// Test two peaks for m/z equality, within a tolerance
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(Peak x, Peak y)
        {
            return _equalityComparer.Equals(x, y);
        }

        /// <inheritdoc />
        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly IEqualityComparer<Peak> _equalityComparer;
    }

    /// <summary>
    /// Compare by m/z; two peaks within ppmTolerance are considered to be equal
    /// </summary>
    public class MzComparerWithPpmTolerance : IEqualityComparer<Peak>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ppmTolerance"></param>
        public MzComparerWithPpmTolerance(double ppmTolerance)
        {
            _ppmTolerance = ppmTolerance;
        }

        /// <summary>
        /// Test two peaks for m/z equality, within a tolerance
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(Peak x, Peak y)
        {
            if (x == null || y == null)
            {
                return false;
            }

            return Math.Abs((x.Mz - y.Mz) / x.Mz * 1E6) <= _ppmTolerance;
        }

        int IEqualityComparer<Peak>.GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly double _ppmTolerance;
    }

    /// <summary>
    /// Compare by m/z; two peaks within toleranceTh Th are considered to be equal
    /// </summary>
    public class MzComparerWithToleranceMz : IEqualityComparer<Peak>
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="toleranceTh"></param>
        public MzComparerWithToleranceMz(double toleranceTh)
        {
            _toleranceTh = toleranceTh;
        }

        /// <summary>
        /// Test two peaks for m/z equality, within a tolerance
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(Peak x, Peak y)
        {
            if (x == null || y == null)
            {
                return false;
            }

            return Math.Abs(x.Mz - y.Mz) <= _toleranceTh;
        }

        /// <inheritdoc />
        public int GetHashCode(Peak p)
        {
            return p.Mz.GetHashCode();
        }

        private readonly double _toleranceTh;
    }

    /// <summary>
    /// Compare by m/z; two peaks within ppmTolerance are considered to be equal
    /// </summary>
    public class MzComparerWithBinning : IComparer<Peak>, IEqualityComparer<Peak>
    {
        /// <summary>
        /// Constructor: 27 bits: max error = 16 ppm, 28 bits (8 ppm), 26 bits (32 ppm)
        /// </summary>
        /// <param name="numBits"></param>
        public MzComparerWithBinning(int numBits = 27)
        {
            NumBits = numBits;
            _numShifts = sizeof(double) * 8 - numBits;
        }

        /// <summary>
        /// Check if two peaks are in the same m/z bin
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(Peak x, Peak y)
        {
            if (x == null || y == null)
            {
                return false;
            }

            return GetBinNumber(x.Mz) == GetBinNumber(y.Mz);
        }

        /// <summary>
        /// Check if two masses are in the same m/z bin
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>True if the items match</returns>
        public bool Equals(double x, double y)
        {
            return GetBinNumber(x) == GetBinNumber(y);
        }

        /// <inheritdoc />
        public int GetHashCode(Peak p)
        {
            return GetBinNumber(p.Mz);
        }

        /// <summary>
        /// Compare two peaks by m/z bin
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns>0, 1, or -1</returns>
        public int Compare(Peak x, Peak y)
        {
            if (x == null)
            {
                return y == null ? 0 : -1;
            }

            if (y == null)
            {
                return 1;
            }

            return GetBinNumber(x.Mz).CompareTo(GetBinNumber(y.Mz));
        }

        /// <summary>
        /// Get the m/z rounded by <see cref="NumBits"/> using bit shifting magic
        /// </summary>
        /// <param name="mz"></param>
        /// <returns>Rounded m/z</returns>
        public double GetRoundedValue(double mz)
        {
            var converted = BitConverter.DoubleToInt64Bits(mz);
            var rounded = (converted >> _numShifts) << _numShifts;
            return BitConverter.Int64BitsToDouble(rounded);
        }

        /// <summary>
        /// Get the bin number for the m/z
        /// </summary>
        /// <param name="mz"></param>
        /// <returns>Bin number</returns>
        public int GetBinNumber(double mz)
        {
            var converted = BitConverter.DoubleToInt64Bits(mz);
            return (int)(converted >> _numShifts);
            //var rounded = (converted >> _numShifts) << _numShifts;
            //return (int)(rounded >> _numShifts);
        }

        /// <summary>
        /// Get the lowest m/z value in bin <paramref name="binNum"/>, inclusive
        /// </summary>
        /// <param name="binNum"></param>
        /// <returns>Lowest m/z in the bin</returns>
        public double GetMzStart(int binNum)
        {
            var rounded = (long)binNum << _numShifts;
            return BitConverter.Int64BitsToDouble(rounded);
        }

        /// <summary>
        /// Get the highest m/z value in bin <paramref name="binNum"/>, exclusive
        /// </summary>
        /// <param name="binNum"></param>
        /// <returns>Highest m/z in the bin</returns>
        public double GetMzEnd(int binNum)
        {
            var rounded = (long)(binNum + 1) << _numShifts;
            return BitConverter.Int64BitsToDouble(rounded);
        }

        /// <summary>
        /// Get average m/z for bin <paramref name="binNum"/>
        /// </summary>
        /// <param name="binNum"></param>
        /// <returns>Average m/z value (halfway between the lowest and highest values)</returns>
        public double GetMzAverage(int binNum)
        {
            return (GetMzStart(binNum) + GetMzEnd(binNum)) * 0.5;
        }

        /// <summary>
        /// Get the ppm tolerance that corresponds to <see cref="NumBits"/>
        /// </summary>
        /// <returns>Tolerance object</returns>
        public Tolerance GetTolerance()
        {
            // 28 bits: max error =  8 ppm
            // 27 bits: max error = 16 ppm
            // 26 bits: max error = 32 ppm

            var numBits = sizeof(double) * 8 - _numShifts;
            var ppm = 16 * Math.Pow(2, 27 - numBits);

            return new Tolerance(ppm);
        }

        /// <summary>
        /// Get the ppm tolerance that corresponds to <see cref="NumBits"/>
        /// </summary>
        public int Ppm
        {
            get
            {
                var numBits = sizeof(double) * 8 - _numShifts;
                var ppm = (int)(16 * Math.Pow(2, 27 - numBits));
                return ppm;
            }
        }

        /// <summary>
        /// Number of bits to use for binning
        /// </summary>
        public readonly int NumBits;
        private readonly int _numShifts;
    }
}
