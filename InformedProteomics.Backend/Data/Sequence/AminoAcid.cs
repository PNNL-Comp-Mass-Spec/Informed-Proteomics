using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class AminoAcid : IMolecule
    {
        #region Constructors

        public AminoAcid(char residue, string name, Composition.Composition comp)
        {
            Residue = residue;
            Name = name;
            Composition = comp;
            Mass = Composition.Mass;
            _nominalMass = Composition.NominalMass;
        }

        #endregion

        #region Properties

        public char Residue { get; private set; }
        public string Name { get; private set; }
        public Composition.Composition Composition { get; private set; }
        public double Mass { get; private set; }

        #endregion

        #region Getters

        //public IsotopomerEnvelope GetIsotopomerEnvelopeRelativeIntensities()
        //{
        //    return Composition.GetIsotopomerEnvelopeRelativeIntensities();
        //}

        public int GetNominalMass()
        {
            return _nominalMass;
        }

        #endregion

        #region Private Members

        private readonly int _nominalMass;

        #endregion

        //public static AminoAcid GetStandardAminoAcid(char residue)
        //{
        //    return StandardAaSet.GetAminoAcid(residue);
        //}

        # region Public Static Members

        public static readonly AminoAcid ProteinNTerm = new AminoAcid('[', "Protein-N-terminus", Data.Composition.Composition.Zero);
        public static readonly AminoAcid ProteinCTerm = new AminoAcid(']', "Protein-C-terminus", Data.Composition.Composition.Zero);

        public static readonly AminoAcid PeptideNTerm = new AminoAcid('(', "Peptide-N-terminus", Data.Composition.Composition.Zero);
        public static readonly AminoAcid PeptideCTerm = new AminoAcid(')', "Peptide-C-terminus", Data.Composition.Composition.Zero);

        public static readonly AminoAcid Empty = new AminoAcid('\0', "Empty", Data.Composition.Composition.Zero);

        public static readonly AminoAcid[] StandardAminoAcidArr =
        {
            new AminoAcid('G', "Glycine", new Composition.Composition(2, 3, 1, 1, 0)),
            new AminoAcid('A', "Alanine", new Composition.Composition(3, 5, 1, 1, 0)),
            new AminoAcid('S', "Serine", new Composition.Composition(3, 5, 1, 2, 0)),
            new AminoAcid('P', "Proline", new Composition.Composition(5, 7, 1, 1, 0)),
            new AminoAcid('V', "Valine", new Composition.Composition(5, 9, 1, 1, 0)),
            new AminoAcid('T', "Threonine", new Composition.Composition(4, 7, 1, 2, 0)),
            new AminoAcid('C', "Cystine", new Composition.Composition(3, 5, 1, 1, 1)),
            new AminoAcid('L', "Leucine", new Composition.Composition(6, 11, 1, 1, 0)),
            new AminoAcid('I', "Isoleucine", new Composition.Composition(6, 11, 1, 1, 0)),
            new AminoAcid('N', "Asparagine", new Composition.Composition(4, 6, 2, 2, 0)),
            new AminoAcid('D', "Aspartate", new Composition.Composition(4, 5, 1, 3, 0)),
            new AminoAcid('Q', "Glutamine", new Composition.Composition(5, 8, 2, 2, 0)),
            new AminoAcid('K', "Lysine", new Composition.Composition(6, 12, 2, 1, 0)),
            new AminoAcid('E', "Glutamate", new Composition.Composition(5, 7, 1, 3, 0)),
            new AminoAcid('M', "Methionine", new Composition.Composition(5, 9, 1, 1, 1)),
            new AminoAcid('H', "Histidine", new Composition.Composition(6, 7, 3, 1, 0)),
            new AminoAcid('F', "Phenylalanine", new Composition.Composition(9, 9, 1, 1, 0)),
            new AminoAcid('R', "Arginine", new Composition.Composition(6, 12, 4, 1, 0)),
            new AminoAcid('Y', "Tyrosine", new Composition.Composition(9, 9, 1, 2, 0)),
            new AminoAcid('W', "Tryptophan", new Composition.Composition(11, 10, 2, 1, 0))
        };

        public const string StandardAminoAcidCharacters = "ACDEFGHIKLMNPQRSTVWY";

        public static bool IsStandardAminoAcidResidue(char residue)
        {
            return StandardAminoAcidCharacters.IndexOf(residue) >= 0;
        }

        //Ala (A) 8.26   Gln (Q) 3.93   Leu (L) 9.66   Ser (S) 6.58
        //Arg (R) 5.53   Glu (E) 6.74   Lys (K) 5.83   Thr (T) 5.34
        //Asn (N) 4.06   Gly (G) 7.08   Met (M) 2.41   Trp (W) 1.09
        //Asp (D) 5.46   His (H) 2.27   Phe (F) 3.86   Tyr (Y) 2.92
        //Cys (C) 1.37   Ile (I) 5.94   Pro (P) 4.71   Val (V) 6.87
        public static readonly double[] StandardAminoAcidFrequency = new double[20]
        {
            0.082674407,0.013712341,0.054649184,0.067460715,0.038634771,
            0.070863777,0.022720448,0.059453508,0.058352517,0.096687018,
            0.02412171,0.040636573,0.047142428,0.039335402,0.055349815,
            0.065859273,0.053448103,0.068761886,0.010909819,0.029226304
        };

        public static double GetUniProtFrequency(char residue)
        {
            var index = StandardAminoAcidCharacters.IndexOf(residue);
            if (index >= 0) return StandardAminoAcidFrequency[index];
            return 0;
        }

        #endregion

        # region Private Static Members

//        private static readonly AminoAcidSet StandardAaSet = new AminoAcidSet();

        #endregion

        public override string ToString()
        {
            return this.Residue.ToString();
        }
    }
}
