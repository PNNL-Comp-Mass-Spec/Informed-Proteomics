using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Amino acid
    /// </summary>
    public class AminoAcid : IMolecule
    {
        // Ignore Spelling: Arg, Asn, Cys, Gln, Glu, Gly, Ile, Leu, Lys, Phe, Ser, Thr, Trp, Tyr, UniProt

        #region Constructors

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="name"></param>
        /// <param name="comp"></param>
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

        /// <summary>
        /// Residue character/symbol
        /// </summary>
        public char Residue { get; }

        /// <summary>
        /// Name of amino acid
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Elemental composition of amino acid
        /// </summary>
        public Composition.Composition Composition { get; }

        /// <summary>
        /// Monoisotopic mass of amino acid
        /// </summary>
        public double Mass { get; }

        #endregion

        #region Getters

        //public IsotopomerEnvelope GetIsotopomerEnvelopeRelativeIntensities()
        //{
        //    return Composition.GetIsotopomerEnvelopeRelativeIntensities();
        //}

        /// <summary>
        /// Nominal mass of the amino acid
        /// </summary>
        /// <returns>Nominal mass</returns>
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

        /// <summary>
        /// Protein N-Terminus marker
        /// </summary>
        public static readonly AminoAcid ProteinNTerm = new('[', "Protein-N-terminus", Data.Composition.Composition.Zero);

        /// <summary>
        /// Protein C-Terminus marker
        /// </summary>
        public static readonly AminoAcid ProteinCTerm = new(']', "Protein-C-terminus", Data.Composition.Composition.Zero);

        /// <summary>
        /// Peptide N-Terminus marker
        /// </summary>
        public static readonly AminoAcid PeptideNTerm = new('(', "Peptide-N-terminus", Data.Composition.Composition.Zero);

        /// <summary>
        /// Peptide C-Terminus marker
        /// </summary>
        public static readonly AminoAcid PeptideCTerm = new(')', "Peptide-C-terminus", Data.Composition.Composition.Zero);

        /// <summary>
        /// Empty amino acid
        /// </summary>
        public static readonly AminoAcid Empty = new('\0', "Empty", Data.Composition.Composition.Zero);

        /// <summary>
        /// Standard amino acids
        /// </summary>
        public static readonly AminoAcid[] StandardAminoAcidArr =
        {
            new AminoAcid('G', "Glycine", new Composition.Composition(2, 3, 1, 1, 0)),
            new AminoAcid('A', "Alanine", new Composition.Composition(3, 5, 1, 1, 0)),
            new AminoAcid('S', "Serine", new Composition.Composition(3, 5, 1, 2, 0)),
            new AminoAcid('P', "Proline", new Composition.Composition(5, 7, 1, 1, 0)),
            new AminoAcid('V', "Valine", new Composition.Composition(5, 9, 1, 1, 0)),
            new AminoAcid('T', "Threonine", new Composition.Composition(4, 7, 1, 2, 0)),
            new AminoAcid('C', "Cysteine", new Composition.Composition(3, 5, 1, 1, 1)),
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

        /// <summary>
        /// String of standard amino acid character symbols
        /// </summary>
        public const string StandardAminoAcidCharacters = "ACDEFGHIKLMNPQRSTVWY";

        /// <summary>
        /// Checks if <paramref name="residue"/> is a standard amino acid character
        /// </summary>
        /// <param name="residue"></param>
        /// <returns>True if a standard amino acid symbol</returns>
        public static bool IsStandardAminoAcidResidue(char residue)
        {
            return StandardAminoAcidCharacters.IndexOf(residue) >= 0;
        }

        /// <summary>
        /// Statistical frequency (per UniProt) of the standard amino acids
        /// </summary>
        public static readonly double[] StandardAminoAcidFrequency = {
            //Ala (A) 8.26   Gln (Q) 3.93   Leu (L) 9.66   Ser (S) 6.58
            //Arg (R) 5.53   Glu (E) 6.74   Lys (K) 5.83   Thr (T) 5.34
            //Asn (N) 4.06   Gly (G) 7.08   Met (M) 2.41   Trp (W) 1.09
            //Asp (D) 5.46   His (H) 2.27   Phe (F) 3.86   Tyr (Y) 2.92
            //Cys (C) 1.37   Ile (I) 5.94   Pro (P) 4.71   Val (V) 6.87
            0.082674407, 0.013712341, 0.054649184, 0.067460715, 0.038634771,
            0.070863777, 0.022720448, 0.059453508, 0.058352517, 0.096687018,
            0.02412171,  0.040636573, 0.047142428, 0.039335402, 0.055349815,
            0.065859273, 0.053448103, 0.068761886, 0.010909819, 0.029226304
        };

        /// <summary>
        /// Get the UniProt statistical frequency of <paramref name="residue"/>
        /// </summary>
        /// <param name="residue"></param>
        /// <returns>Frequency at which this residue is seen (number between 0 and 1)</returns>
        public static double GetUniProtFrequency(char residue)
        {
            var index = StandardAminoAcidCharacters.IndexOf(residue);
            if (index >= 0)
            {
                return StandardAminoAcidFrequency[index];
            }

            return 0;
        }

        #endregion

        /// <inheritdoc />
        public override string ToString()
        {
            return Residue.ToString();
        }
    }
}
