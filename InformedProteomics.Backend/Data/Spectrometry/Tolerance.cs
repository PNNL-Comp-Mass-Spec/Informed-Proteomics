using System;
using System.Linq;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Facilitates working with mass spec measurement tolerances
    /// </summary>
    public class Tolerance
    {
        // Ignore Spelling: Da, Daltons, Thomsons

        private readonly double _value;
        private readonly ToleranceUnit _unit;

        /// <summary>
        /// Constant tolerance object for 1 ppm
        /// </summary>
        public static readonly Tolerance OnePpm = new(1);

        /// <summary>
        /// Instantiate a tolerance with the supplied value and ppm units
        /// </summary>
        /// <param name="value"></param>
        public Tolerance(double value)
            : this(value, ToleranceUnit.Ppm)
        {
        }

        /// <summary>
        /// Instantiate a tolerance with the supplied value and units
        /// </summary>
        /// <param name="value"></param>
        /// <param name="unit"></param>
        public Tolerance(double value, ToleranceUnit unit)
        {
            _value = value;
            _unit = unit;
        }

        /// <summary>
        /// Get the tolerance value for this instance
        /// </summary>
        /// <returns>Tolerance value</returns>
        public double GetValue()
        {
            return _value;
        }

        /// <summary>
        /// Get the tolerance unit set for this instance
        /// </summary>
        /// <returns>Tolerance unit</returns>
        public ToleranceUnit GetUnit()
        {
            return _unit;
        }

        /// <summary>
        /// Get the tolerance limit in terms of m/z for the supplied mass
        /// </summary>
        /// <param name="mz"></param>
        /// <returns>Tolerance, as m/z</returns>
        public double GetToleranceAsMz(double mz)
        {
            if (_unit == ToleranceUnit.Mz)
            {
                return _value;
            }

            if (_unit == ToleranceUnit.Da)
            {
                throw new System.ArgumentException("This function cannot be called for Da unit");
            }

            return mz * _value / 1E6;
        }

        /// <summary>
        /// Get the tolerance limit in terms of Thomsons (m/z) for the supplied mass
        /// </summary>
        /// <param name="mz"></param>
        /// <returns>Tolerance, as m/z</returns>
        [Obsolete("Use GetToleranceAsMz")]
        public double GetToleranceAsTh(double mz)
        {
            return GetToleranceAsMz(mz);
        }

        /// <summary>
        /// Get the tolerance limit in terms of Daltons for the supplied mass and charge
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="charge"></param>
        /// <returns>Tolerance, as uncharged mass in Daltons</returns>
        public double GetToleranceAsDa(double mz, int charge)
        {
            if (_unit == ToleranceUnit.Da)
            {
                return _value;
            }

            return GetToleranceAsMz(mz) * charge;
        }

        /// <summary>
        /// Test the supplied m/z values to see if they are within the tolerance limits in this instance
        /// </summary>
        /// <param name="mz1"></param>
        /// <param name="mz2"></param>
        /// <returns>True if the m/z values are within tolerance</returns>
        public bool IsWithin(double mz1, double mz2)
        {
            var tolTh = GetToleranceAsMz(mz1);
            if (tolTh < 0)
            {
                tolTh = -tolTh;
            }

            return mz2 > mz1 - tolTh && mz2 < mz1 + tolTh;
        }

        /// <summary>
        /// Return the string representation of this tolerance
        /// </summary>
        /// <returns>Human-readable tolerance, including units</returns>
        public override string ToString()
        {
            return string.Format("{0}{1}", _value, _unit);
        }

        /// <summary>
        /// Return a tolerance object that was created using the data from the supplied string
        /// </summary>
        /// <param name="tolStr"></param>
        /// <returns>Tolerance object</returns>
        public static Tolerance Parse(string tolStr)
        {
            tolStr = tolStr.ToLower();
            var units = System.Enum.GetValues(typeof(ToleranceUnit)).Cast<ToleranceUnit>();
            var unit = units.FirstOrDefault(u => tolStr.Contains(u.ToString().ToLower()));
            var index = tolStr.IndexOf(unit.ToString().ToLower(), StringComparison.Ordinal);
            if (index < 0)
            {
                return null;
            }

            var valueStr = tolStr.Substring(0, index);
            var unitStr = tolStr.Substring(index, tolStr.Length - index);
            if (string.IsNullOrEmpty(valueStr) || string.IsNullOrEmpty(unitStr))
            {
                return null;
            }

            ToleranceUnit tolUnit;
            switch (unitStr)
            {
                case "da":
                    tolUnit = ToleranceUnit.Da;
                    break;
                case "th":
                case "mz":
                case "m/z":
                    tolUnit = ToleranceUnit.Mz;
                    break;
                case "ppm":
                default:
                    tolUnit = ToleranceUnit.Ppm;
                    break;
            }

            return new Tolerance(Convert.ToDouble(valueStr), tolUnit);
        }
    }
}
