using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Constants
    {
        public const double RescalingConstant = 0.9995;
        public const double RescalingConstantHighPrecision = 274.335215;
        //public static readonly double H2O = Annotation.Composition.H2O.GetMass();
        public const double C13 = 13.00335483;
        public const double C13MinusC12 = C13 - 12.0;
        public const double Proton = 1.00727649;

        public static int GetBinNum(double m)
        {
            return (int)Math.Round(m * RescalingConstant);
        }

        public static int GetBinNumHighPrecision(double m)
        {
            return (int) Math.Round(m*RescalingConstantHighPrecision);
        }
    }
}
