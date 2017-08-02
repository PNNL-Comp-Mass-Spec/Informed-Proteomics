using System;
using System.Collections.Generic;
using System.Diagnostics;
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
            Composition = modificationList.Aggregate(Data.Composition.Composition.Zero, (current, mod) => current + mod.Composition);
        }

        public Composition.Composition Composition { get; }
        public IList<Modification> Modifications { get; }

        public int GetNumModifications()
        {
            return Modifications.Count;
        }

        public override string ToString()
        {
            Debug.Assert(Modifications != null, "Modifications != null");
            return string.Join(",",Modifications);
        }
    }
}
