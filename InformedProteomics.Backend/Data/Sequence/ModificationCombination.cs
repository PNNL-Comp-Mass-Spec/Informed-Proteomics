using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Combination of modifications applied to an amino acid sequence
    /// </summary>
    public class ModificationCombination
    {
        /// <summary>
        /// No modification
        /// </summary>
        public static readonly ModificationCombination NoModification = new ModificationCombination(new List<Modification>());

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="modificationList"></param>
        public ModificationCombination(IList<Modification> modificationList)
        {
            Modifications = modificationList;
            Composition = modificationList.Aggregate(Data.Composition.Composition.Zero, (current, mod) => current + mod.Composition);
        }

        /// <summary>
        /// Composition
        /// </summary>
        public Composition.Composition Composition { get; }

        /// <summary>
        /// Modifications applied to the sequence
        /// </summary>
        public IList<Modification> Modifications { get; }

        /// <summary>
        /// Number of modifications applied to the sequence
        /// </summary>
        /// <returns>Modification count</returns>
        public int GetNumModifications()
        {
            return Modifications.Count;
        }

        /// <inheritdoc />
        public override string ToString()
        {
            Debug.Assert(Modifications != null, "Modifications != null");
            return string.Join(",", Modifications);
        }
    }
}
