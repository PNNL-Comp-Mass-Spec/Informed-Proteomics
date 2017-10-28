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
        /// <param name="suffixComposition"></param>
        public Cleavage(Composition.Composition prefixComposition, Composition.Composition suffixComposition)
        {
            PrefixComposition = prefixComposition;
            SuffixComposition = suffixComposition;
        }

        /// <summary>
        /// Prefix composition
        /// </summary>
        public Composition.Composition PrefixComposition { get; }

        /// <summary>
        /// Suffix composition
        /// </summary>
        public Composition.Composition SuffixComposition { get; }
    }
}
