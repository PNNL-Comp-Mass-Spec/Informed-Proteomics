using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Biology
{
    class Constants
    {
        public static readonly double RescalingConstant = 0.9995;
        public static readonly double H2O = Sequence.Composition.H2O.GetMass();
        public static readonly double H = Atom.Get("H").Mass;
    }
}
