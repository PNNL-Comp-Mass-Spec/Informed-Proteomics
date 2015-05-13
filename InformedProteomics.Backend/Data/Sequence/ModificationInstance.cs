namespace InformedProteomics.Backend.Data.Sequence
{
    public class ModificationInstance
    {
        public ModificationInstance(Modification modification, int index)
        {
            Modification = modification;
            Index = index;
        }

        public Modification Modification { get; private set; }
        public int Index { get; private set; }

        public override string ToString()
        {
            return Modification.Name + " " + (Index+1);
        }

        public ModificationInstance GetModificationInstanceWithOffset(int offset)
        {
            return new ModificationInstance(Modification, Index + offset);
        }
    }
}
