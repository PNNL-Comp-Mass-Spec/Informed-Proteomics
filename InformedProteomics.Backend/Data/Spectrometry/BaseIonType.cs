using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class BaseIonType
    {
        private BaseIonType(string symbol, bool isPrefix, Composition offsetComposition)
        {
            Symbol = symbol;
            IsPrefix = isPrefix;
            OffsetComposition = offsetComposition;
        }

        public string Symbol { get; private set; }
        public bool IsPrefix { get; private set; }
        public Composition OffsetComposition { get; private set; }

        public static readonly BaseIonType A = new BaseIonType("a", true, new Composition(-1, 0, 0, -1, 0)); // -CO
        public static readonly BaseIonType B = new BaseIonType("b", true, Composition.Zero);  
        public static readonly BaseIonType C = new BaseIonType("c", true, Composition.NH3); 

        public static readonly BaseIonType X = new BaseIonType("x", false, Composition.H2O+Composition.CO);
        public static readonly BaseIonType Y = new BaseIonType("y", false, Composition.H2O);
        public static readonly BaseIonType Z = new BaseIonType("x", false, Composition.H2O-Composition.NH2);

        public static readonly IEnumerable<BaseIonType> AllBaseIonTypes = new List<BaseIonType> {A, B, C, X, Y, Z};
    }
}
