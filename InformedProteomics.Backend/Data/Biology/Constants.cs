using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Constants
    {
        public static readonly double RescalingConstant = 0.9995;
        public static readonly double H2O = Sequence.Composition.H2O.GetMass();
        public static readonly double C13 = 13.00335483;
        public static readonly double C13MinusC12 = C13 - 12.0;
        public static readonly double Proton = 1.00727649;
    }
}
