namespace InformedProteomics.Backend.Data.Spectrometry
{
    using System;
    using System.Collections.Generic;
    using System.Linq;

    public class Tolerance
    {
        private readonly double _value;
        private readonly ToleranceUnit _unit;

        public static readonly Tolerance OnePpm = new Tolerance(1);

        public Tolerance(double value)
            : this(value, ToleranceUnit.Ppm)
        {
        }

        public Tolerance(double value, ToleranceUnit unit)
        {
            _value = value;
            _unit = unit;
        }

        public double GetValue()
        {
            return _value;
        }

        public ToleranceUnit GetUnit()
        {
            return _unit;
        }

        public double GetToleranceAsTh(double mz)
        {
            if (_unit == ToleranceUnit.Th) return _value;
            if (_unit == ToleranceUnit.Da)
                throw new System.ArgumentException("This function cannot be called for Da unit");
            return mz * _value / 1E6;
        }

        public double GetToleranceAsDa(double mz, int charge)
        {
            if (_unit == ToleranceUnit.Da) return _value;
            return GetToleranceAsTh(mz)*charge;
        }

        public bool IsWithin(double mz1, double mz2)
        {
            var tolTh = GetToleranceAsTh(mz1);
            if (tolTh < 0) tolTh = -tolTh;
            return mz2 > mz1 - tolTh && mz2 < mz1 + tolTh;
        }

        public override string ToString()
        {
            return string.Format("{0}{1}", _value, _unit);
        }

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

            return new Tolerance(Convert.ToDouble(valueStr), (ToleranceUnit)Enum.Parse(typeof(ToleranceUnit), unitStr));
        }
    }
}
