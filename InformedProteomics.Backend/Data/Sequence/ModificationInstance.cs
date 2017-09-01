namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Instance of a modification
    /// </summary>
    public class ModificationInstance
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="modification"></param>
        /// <param name="index"></param>
        public ModificationInstance(Modification modification, int index)
        {
            Modification = modification;
            Index = index;
        }

        /// <summary>
        /// Modification
        /// </summary>
        public Modification Modification { get; }

        /// <summary>
        /// Index of the modification
        /// </summary>
        public int Index { get; }

        /// <inheritdoc />
        public override string ToString()
        {
            return Modification.Name + " " + (Index+1);
        }

        /// <summary>
        /// Get a new <see cref="ModificationInstance"/> at offset <paramref name="offset"/> from this instance
        /// </summary>
        /// <param name="offset"></param>
        /// <returns></returns>
        public ModificationInstance GetModificationInstanceWithOffset(int offset)
        {
            return new ModificationInstance(Modification, Index + offset);
        }
    }
}
