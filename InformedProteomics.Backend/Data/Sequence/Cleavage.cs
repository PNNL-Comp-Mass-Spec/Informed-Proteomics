namespace InformedProteomics.Backend.Data.Sequence
{
    public class Cleavage
    {
        public Cleavage(Composition.Composition prefixComposition, Composition.Composition suffixComposition)
        {
            PrefixComposition = prefixComposition;
            SuffixComposition = suffixComposition;
        }

        public Composition.Composition PrefixComposition { get; private set; }
        public Composition.Composition SuffixComposition { get; private set; }
    }
}
