using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class NeutralLoss
    {
        public NeutralLoss(string name, string symbol, Composition.Composition composition)
        {
            Symbol = symbol;
            Name = name;
            Composition = composition;
        }

        public string Symbol { get; private set; }
        public string Name { get; private set; }
        public Composition.Composition Composition { get; private set; }

        public static readonly NeutralLoss NoLoss = new NeutralLoss("", "NoLoss", Data.Composition.Composition.Zero);
        public static readonly NeutralLoss H2O = new NeutralLoss("-H2O", "H2O", Data.Composition.Composition.H2O);
        public static readonly NeutralLoss NH3 = new NeutralLoss("-NH3", "NH3", Data.Composition.Composition.NH3);

        //public static readonly NeutralLoss DeconvolutedIon = new NeutralLoss("'", new CompositionWithDeltaMass(Biology.Constants.Proton));

        public static readonly IEnumerable<NeutralLoss> CommonNeutralLosses = new List<NeutralLoss> { NoLoss, H2O, NH3 };
    }
}
