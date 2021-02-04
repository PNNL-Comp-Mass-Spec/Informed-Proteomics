using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// A modification specified in the search
    /// </summary>
    public class SearchModification
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mod"></param>
        /// <param name="targetResidue"></param>
        /// <param name="loc"></param>
        /// <param name="isFixedModification"></param>
        public SearchModification(Modification mod, char targetResidue, SequenceLocation loc, bool isFixedModification)
        {
            Modification = mod;
            TargetResidue = targetResidue;
            Location = loc;
            IsFixedModification = isFixedModification;
        }

        /// <summary>
        /// Empty constructor
        /// </summary>
        public SearchModification()
        {
        }

        /// <summary>
        /// Modification
        /// </summary>
        public Modification Modification { get; set; }

        /// <summary>
        /// Residue targeted by <see cref="Modification"/>
        /// </summary>
        public char TargetResidue { get; set; }

        /// <summary>
        /// Allowed location in sequence of <see cref="Modification"/>
        /// </summary>
        public SequenceLocation Location { get; set; }

        /// <summary>
        /// If the modification is fixed/static
        /// </summary>
        public bool IsFixedModification { get; set; }

        /// <summary>
        /// Name of the modification
        /// </summary>
        public string Name => Modification.Name;

        /// <summary>
        /// Mass of the modification
        /// </summary>
        public double Mass => Modification.Mass;

        /// <inheritdoc />
        public override string ToString()
        {
            return ToString(false);
        }

        /// <summary>
        /// String representation of this search mod
        /// </summary>
        /// <param name="spaceSeparated">When true, separate fields using spaces instead of commas</param>
        /// <returns></returns>
        public string ToString(bool spaceSeparated)
        {
            var formatString = spaceSeparated ?
                                   "{0,-30} {1,-5} {2,-4} {3,-15} {4}" :
                                   "{0},{1},{2},{3},{4}";

            return string.Format(formatString,
                Modification.Composition.ToString(spaceSeparated),
                TargetResidue,
                IsFixedModification ? "fix" : "opt",
                Location,
                Modification.Name
            );
        }
    }
}
