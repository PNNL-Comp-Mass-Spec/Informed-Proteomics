using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UIMFLibrary;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Tolerance
    {
        private readonly double _value;
        private readonly DataReader.ToleranceType _unit;

        public Tolerance(double value) : this(value, DataReader.ToleranceType.PPM)
        {
        }

        public Tolerance(double value, DataReader.ToleranceType unit)
        {
            _value = value;
            _unit = unit;
        }

        public double GetValue()
        {
            return _value;
        }

        public DataReader.ToleranceType GetUnit()
        {
            return _unit;
        }

        public double GetToleranceAsTh(double mz)
        {
            if (_unit == DataReader.ToleranceType.Thomson)
            {
                return _value;
            }
            else
            {
                return mz*_value/1E6;
            }
        }
    }
}
