namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Cleavage of a composition
    /// </summary>
    public class Cleavage
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="prefixComposition"></param>
        /// <param name="prefixResidue"></param>
        /// <param name="suffixComposition"></param>
        /// <param name="suffixResidue"></param>
        public Cleavage(Composition.Composition prefixComposition, AminoAcid prefixResidue, Composition.Composition suffixComposition, AminoAcid suffixResidue)
        {
            PrefixComposition = prefixComposition;
            SuffixComposition = suffixComposition;

            PrefixResidue = prefixResidue;
            SuffixResidue = suffixResidue;
        }

        /// <summary>
        /// Prefix composition
        /// </summary>
        public Composition.Composition PrefixComposition { get; }

        /// <summary>
        /// Suffix composition
        /// </summary>
        public Composition.Composition SuffixComposition { get; }

        /// <summary>
        /// Prefix residue
        /// </summary>
        public AminoAcid PrefixResidue { get; }

        /// <summary>
        /// Suffix residue
        /// </summary>
        public AminoAcid SuffixResidue { get; set; }
    }
}
