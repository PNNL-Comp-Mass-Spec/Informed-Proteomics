using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Combination of modifications applied to an amino acid sequence
    /// </summary>
    public class ModificationCombination
    {
        public static readonly ModificationCombination NoModification = new ModificationCombination(new List<Modification>());

        public ModificationCombination(IList<Modification> modificationList)
        {
            Modifications = modificationList;
            Composition = modificationList.Aggregate(Composition.Zero, (current, mod) => current + mod.Composition);
        }

        public Composition Composition { get; private set; }
        public IList<Modification> Modifications { get; private set; }
    }
}
