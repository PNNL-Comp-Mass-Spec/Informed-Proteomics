using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class BaseIonType
    {
        private BaseIonType(string symbol, bool isPrefix, Composition.Composition offsetComposition)
        {
            Symbol = symbol;
            IsPrefix = isPrefix;
            OffsetComposition = offsetComposition;
        }

        public string Symbol { get; private set; }
        public bool IsPrefix { get; private set; }
        public Composition.Composition OffsetComposition { get; private set; }

        public static readonly BaseIonType A = new BaseIonType("a", true, new Composition.Composition(-1, 0, 0, -1, 0)); // -CO
        public static readonly BaseIonType B = new BaseIonType("b", true, Composition.Composition.Zero);  
        public static readonly BaseIonType C = new BaseIonType("c", true, Composition.Composition.NH3); 

        public static readonly BaseIonType X = new BaseIonType("x", false, Composition.Composition.H2O+Composition.Composition.CO);
        public static readonly BaseIonType Y = new BaseIonType("y", false, Composition.Composition.H2O);

        // Z. instead of Z
        public static readonly BaseIonType Z = new BaseIonType("z", false, Composition.Composition.H2O-Composition.Composition.NH2+Composition.Composition.Hydrogen);  

        public static readonly IEnumerable<BaseIonType> AllBaseIonTypes = new List<BaseIonType> {A, B, C, X, Y, Z};
    }
}
