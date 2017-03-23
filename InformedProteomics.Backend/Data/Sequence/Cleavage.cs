namespace InformedProteomics.Backend.Data.Sequence
{
    public class Cleavage
    {
        public Cleavage(Composition.Composition prefixComposition, AminoAcid prefixResidue, Composition.Composition suffixComposition, AminoAcid suffixResidue)
        {
            PrefixComposition = prefixComposition;
            SuffixComposition = suffixComposition;

            PrefixResidue = prefixResidue;
            SuffixResidue = suffixResidue;
        }

        public Composition.Composition PrefixComposition { get; private set; }
        public Composition.Composition SuffixComposition { get; private set; }

        public AminoAcid PrefixResidue { get; private set; }

        public AminoAcid SuffixResidue { get; set; }
    }
}
