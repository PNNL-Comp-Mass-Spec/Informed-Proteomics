using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class SearchModification
    {
        public SearchModification(Modification mod, AminoAcid targetAA, SequenceLocation loc, bool isFixedModification)
        {
            Modification = mod;
            TargetAA = targetAA;
            Location = loc;
            IsFixedModification = isFixedModification;
        }

        public Modification Modification { get; private set; }
        public AminoAcid TargetAA { get; private set; }
        public SequenceLocation Location { get; private set; }
        public bool IsFixedModification { get; private set; }
    }
}
