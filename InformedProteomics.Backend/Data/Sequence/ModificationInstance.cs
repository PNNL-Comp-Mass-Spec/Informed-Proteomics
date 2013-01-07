using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class ModificationInstance
    {
        public ModificationInstance(Modification mod, AminoAcid targetAA, SequenceLocation loc)
        {
            Modification = mod;
            TargetAA = targetAA;
            Location = loc;
        }

        public Modification Modification { get; private set; }
        public AminoAcid TargetAA { get; private set; }
        public SequenceLocation Location { get; private set; }
    }
}
