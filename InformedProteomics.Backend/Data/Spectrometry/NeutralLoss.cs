using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Neutral loss
    /// </summary>
    public class NeutralLoss
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name"></param>
        /// <param name="symbol"></param>
        /// <param name="composition"></param>
        public NeutralLoss(string name, string symbol, Composition.Composition composition)
        {
            Symbol = symbol;
            Name = name;
            Composition = composition;
        }

        /// <summary>
        /// Symbol for the type of neutral loss
        /// </summary>
        public string Symbol { get; }

        /// <summary>
        /// Name of the neutral loss type
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Neutral loss composition
        /// </summary>
        public Composition.Composition Composition { get; }

        /// <summary>
        /// No neutral loss
        /// </summary>
        public static readonly NeutralLoss NoLoss = new(string.Empty, "NoLoss", Data.Composition.Composition.Zero);

        /// <summary>
        /// Neutral water loss
        /// </summary>
        public static readonly NeutralLoss H2O = new("-H2O", "H2O", Data.Composition.Composition.H2O);

        /// <summary>
        /// Neutral ammonia loss
        /// </summary>
        public static readonly NeutralLoss NH3 = new("-NH3", "NH3", Data.Composition.Composition.NH3);

        //public static readonly NeutralLoss DeconvolutedIon = new NeutralLoss("'", new CompositionWithDeltaMass(Biology.Constants.Proton));

        /// <summary>
        /// List of common neutral losses
        /// </summary>
        public static readonly IEnumerable<NeutralLoss> CommonNeutralLosses = new List<NeutralLoss> { NoLoss, H2O, NH3 };
    }
}
