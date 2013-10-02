using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UIMFLibrary;

namespace InformedProteomics.Backend.Data.Spectrometry
{
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
            if (_unit == ToleranceUnit.Thomson) return _value;
            if (_unit == ToleranceUnit.Dalton)
                throw new System.ArgumentException("This function cannot be called for Dalton unit");
            return mz * _value / 1E6;
        }

        public double GetToleranceAsDa(double mz, int charge)
        {
            if (_unit == ToleranceUnit.Dalton) return _value;
            return GetToleranceAsTh(mz)*charge;
        }
    }
}
