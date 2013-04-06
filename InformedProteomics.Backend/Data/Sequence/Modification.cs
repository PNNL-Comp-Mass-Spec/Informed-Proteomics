using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Modification : IMolecule
    {
        public int AccessionNum { get; private set; }
        public Composition Composition { get; private set; }
        public string Name { get; private set; }

        public Composition GetComposition()
        {
            return Composition;
        }

        public double GetMass()
        {
            return Composition.GetMass();
        }

        private Modification(int accessionNum, Composition composition, string name)
        {
            AccessionNum = accessionNum;
            Composition = composition;
            Name = name;
        }

        public override int GetHashCode()
        {
            return AccessionNum;
        }

        public override bool Equals(object obj)
        {
            var otherMod = obj as Modification;
            return otherMod != null && AccessionNum == otherMod.AccessionNum;
        }

        /// <summary>
        /// Returns a string that represents this modification object.
        /// </summary>
        /// <returns>
        /// A string that represents the modification.
        /// </returns>
        public override string ToString()
        {
            return Name;
        }

        public static Modification Get(string psiMsName)
        {
            return NameToModMap[psiMsName];
        }

        public static readonly Modification NoModification = new Modification(0, new Composition(0, 0, 0, 0, 0), "No modification");
        public static readonly Modification Carbamidomethylation = new Modification(4, new Composition(2, 3, 1, 1, 0), "Carbamidomethyl");
        public static readonly Modification Phosphorylation = new Modification(21, new Composition(0, 1, 0, 3, 0, 1), "Phospho");
        public static readonly Modification Oxidation = new Modification(35, new Composition(0, 0, 0, 1, 0), "Oxidation");
        public static readonly Modification Acetylation = new Modification(1, new Composition(2, 2, 0, 1, 0), "Acetyl");
        public static readonly Modification Deamidation = new Modification(7, new Composition(0, -1, -1, 1, 0), "Deamidated");
        public static readonly Modification PyroGluQ = new Modification(28, new Composition(0, -3, -1, 0, 0), "Gln->pyro-Glu");

        private static readonly Dictionary<string, Modification> NameToModMap;

        private static readonly Modification[] CommonModifications =
            {
                Acetylation,
                Carbamidomethylation,
                new Modification(5, new Composition(1, 1, 1, 1, 0), "Carbamyl"),
                new Modification(6, new Composition(2, 2, 2, 0, 0), "Carboxymethyl"),
                Deamidation,
                new Modification(17, new Composition(5, 9, 1, 1, 0), "NIPCAM"),
                Phosphorylation,
                new Modification(26, new Composition(0, -3, -1, 0, 0), "Pyro-carbamidomethyl"),
                new Modification(27, new Composition(0, -2, 0, -1, 0), "Glu->pyro-Glu"),
                PyroGluQ,
                new Modification(34, new Composition(1, 2, 0, 0, 0), "Methyl"),
                Oxidation
            };

        static Modification()
        {
            NameToModMap = new Dictionary<string, Modification>();
            foreach (var modification in CommonModifications)
            {
                NameToModMap.Add(modification.Name, modification);
            }
        }
    }
}
