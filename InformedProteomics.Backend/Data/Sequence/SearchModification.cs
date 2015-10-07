using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// A modification specified in the search
    /// </summary>
    public class SearchModification
    {
        public SearchModification(Modification mod, char targetResidue, SequenceLocation loc, bool isFixedModification)
        {
            Modification = mod;
            TargetResidue = targetResidue;
            Location = loc;
            IsFixedModification = isFixedModification;
        }

        public SearchModification()
        {
            
        }

        public Modification Modification { get; set; }
        public char TargetResidue { get; set; }
        public SequenceLocation Location { get; set; }
        public bool IsFixedModification { get; set; }
        public override string ToString()
        {
            return string.Format("{0},{1},{2},{3},{4}", 
                Modification.Composition, 
                TargetResidue, 
                (IsFixedModification ? "fix" : "opt"), 
                Location, 
                Modification.Name
                );
        }
    }
}
