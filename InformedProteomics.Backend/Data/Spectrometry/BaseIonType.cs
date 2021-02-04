using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using Composition = InformedProteomics.Backend.Data.Composition.Composition;

    /// <summary>
    /// Base class for IonTypes
    /// </summary>
    public sealed class BaseIonType
    {
        /// <summary>
        /// Default base ion types
        /// </summary>
        public static readonly BaseIonType A, Ar, B, C, D, V, W, X, Xr, Y, YM1, Zr, Z;

        /// <summary>
        /// List of all base ion types
        /// </summary>
        public static readonly IEnumerable<BaseIonType> AllBaseIonTypes;

        private readonly CompositionCalculator _compositionCalculator;

        private BaseIonType(string symbol, bool isPrefix, Composition offsetComposition, CompositionCalculator compositionCalculator = null)
        {
            Symbol = symbol;
            IsPrefix = isPrefix;
            OffsetComposition = offsetComposition;
            _compositionCalculator = compositionCalculator ?? (_ => new List<Composition> { OffsetComposition });
        }

        static BaseIonType()
        {
            // For A, B, C, X, Y, Z ions: sum of all amino acids + offset
            A = new BaseIonType("a", true, new Composition(-1, 0, 0, -1, 0)); // -CO
            Ar = new BaseIonType("a.", true, A.OffsetComposition + Composition.Hydrogen);
            B = new BaseIonType("b", true, Composition.Zero);
            C = new BaseIonType("c", true, Composition.NH3);

            //X = new BaseIonType("x", false, Composition.H2O + Composition.CO);
            X = new BaseIonType("x", false, new CompositionWithDeltaMass(44.9977) - Composition.Hydrogen);
            Xr = new BaseIonType("x.", false, X.OffsetComposition + Composition.Hydrogen);
            Y = new BaseIonType("y", false, Composition.H2O);
            YM1 = new BaseIonType("y-1", false, Y.OffsetComposition - Composition.Hydrogen);

            Z = new BaseIonType("z", false, Composition.H2O - Composition.NH2);
            Zr = new BaseIonType("z.", false, Z.OffsetComposition + Composition.Hydrogen);

            var aminoAcidSet = new AminoAcidSet();

            // D ions have additional options for isoleucine and threonine. All offsets defined as sum of previous residues + offset.
            var dIonOffsets = new Dictionary<char, double[]>
            {   // Defined in terms of sum of previous residue weights + offset
                { '*', new [] { 44.0500 - Composition.Hydrogen.Mass } },    // For all residues except for V, I and T
                { 'V', new [] { 58.0657 - Composition.Hydrogen.Mass } },
                { 'I', new [] { 58.0657 - Composition.Hydrogen.Mass, 72.0813 - Composition.Hydrogen.Mass } }, // for isoleucine
                { 'T', new [] { 58.0657 - Composition.Hydrogen.Mass, 60.0450 - Composition.Hydrogen.Mass } }  // for threonine
            };
            D = new BaseIonType("d", true, new CompositionWithDeltaMass(44.0500), aminoAcid => DefaultCompositionCalculator(dIonOffsets, aminoAcid));

            // V only has one option for all terminal residues. Sum of previous residues + offset
            V = new BaseIonType("v", false, new CompositionWithDeltaMass(74.0242), r => new List<Composition> { r == null ? new CompositionWithDeltaMass(74.0242) : new CompositionWithDeltaMass(74.0242) - r.Composition - Composition.Hydrogen });

            // W ions have additional options for isoleucine and threonine. All offsets defined as sum of previous residues + offset.
            var wIonOffsets = new Dictionary<char, double[]>
            {   // Defined in terms of sum of previous residue weights + offset
                { '*', new [] { 73.0290 - Composition.Hydrogen.Mass } },    // For all residues except for V, I and T
                { 'V', new [] { 87.0446 - Composition.Hydrogen.Mass } },
                { 'I', new [] { 87.0446 - Composition.Hydrogen.Mass, aminoAcidSet.GetAminoAcid('I').Mass - 12.0238 - Composition.Hydrogen.Mass } }, // for isoleucine
                { 'T', new [] { 87.0446 - Composition.Hydrogen.Mass, aminoAcidSet.GetAminoAcid('T').Mass - 12.0238 - Composition.Hydrogen.Mass } }  // for threonine
            };
            W = new BaseIonType("w", false, new CompositionWithDeltaMass(73.0290), aminoAcid => DefaultCompositionCalculator(wIonOffsets, aminoAcid));

            AllBaseIonTypes = new List<BaseIonType> { A, Ar, B, C, D, V, W, X, Xr, Y, YM1, Z, Zr };
        }

        /// <summary>
        /// Calculates the composition of the ion, with option amino acid
        /// </summary>
        /// <param name="aminoAcid"></param>
        /// <returns></returns>
        public delegate IEnumerable<Composition> CompositionCalculator(AminoAcid aminoAcid = null);

        /// <summary>
        /// Ion symbol
        /// </summary>
        public string Symbol { get; }

        /// <summary>
        /// If the ion is a prefix ion
        /// </summary>
        public bool IsPrefix { get; }

        /// <summary>
        /// The offset this ion adds to a composition
        /// </summary>
        public Composition OffsetComposition { get; }

        /// <summary>
        /// Gets the deconvoluted version of this ion
        /// </summary>
        /// <returns></returns>
        public BaseIonType GetDeconvolutedIon()
        {
            return new BaseIonType(Symbol + "'", IsPrefix,
                OffsetComposition + new CompositionWithDeltaMass(-Biology.Constants.Proton));
        }

        /// <summary>
        /// Get possible compositions of this ion
        /// </summary>
        /// <param name="residue"></param>
        /// <returns></returns>
        public IEnumerable<Composition> GetPossibleCompositions(AminoAcid residue = null)
        {
            return _compositionCalculator(residue);
        }

        private static IEnumerable<Composition> DefaultCompositionCalculator(
            IReadOnlyDictionary<char, double[]> possibleOffsets,
            AminoAcid aminoAcid)
        {
            if (aminoAcid == null)
            {
                return possibleOffsets['*'].Select(off => new CompositionWithDeltaMass(off));
            }

            var key = possibleOffsets.ContainsKey(aminoAcid.Residue) ? aminoAcid.Residue : '*';
            return possibleOffsets[key].Select(offset => new CompositionWithDeltaMass(offset) - aminoAcid.Composition);
        }
    }
}
