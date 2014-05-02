using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Modification : IMolecule
    {
        public int AccessionNum { get; private set; }
        public Composition.Composition Composition { get; private set; }
        public string Name { get; private set; }

        public Composition.Composition GetComposition()
        {
            return Composition;
        }

        public double GetMass()
        {
            return Composition.Mass;
        }

        private Modification(int accessionNum, Composition.Composition composition, string name)
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

        public static IList<Modification> GetFromMass(string mass)
        {
            if (mass.StartsWith("+")) mass = mass.Substring(1);
            IList<Modification> modList;
            return MassToModMap.TryGetValue(mass, out modList) ? modList : null;
        }

        public static readonly Modification NoModification = new Modification(0, new Composition.Composition(0, 0, 0, 0, 0), "No modification");

        public static readonly Modification Acetylation = new Modification(1, new Composition.Composition(2, 2, 0, 1, 0), "Acetyl");
        public static readonly Modification Carbamidomethylation = new Modification(4, new Composition.Composition(2, 3, 1, 1, 0), "Carbamidomethyl");
        public static readonly Modification Carbamylation = new Modification(5, new Composition.Composition(1, 1, 1, 1, 0), "Carbamyl");
        public static readonly Modification Carboxymethylation = new Modification(6, new Composition.Composition(2, 2, 2, 0, 0), "Carboxymethyl");
        public static readonly Modification Deamidation = new Modification(7, new Composition.Composition(0, -1, -1, 1, 0), "Deamidated");
        public static readonly Modification NipCam = new Modification(17, new Composition.Composition(5, 9, 1, 1, 0), "NIPCAM");
        public static readonly Modification Phosphorylation = new Modification(21, new Composition.Composition(0, 1, 0, 3, 0, 1), "Phospho");
        public static readonly Modification PyroCarbamidomethyl = new Modification(26, new Composition.Composition(0, -3, -1, 0, 0), "Pyro-carbamidomethyl");
        public static readonly Modification PyroGluE = new Modification(27, new Composition.Composition(0, -2, 0, -1, 0), "Glu->pyro-Glu");
        public static readonly Modification PyroGluQ = new Modification(28, new Composition.Composition(0, -3, -1, 0, 0), "Gln->pyro-Glu");
        public static readonly Modification Methylation = new Modification(34, new Composition.Composition(1, 2, 0, 0, 0), "Methyl");
        public static readonly Modification Oxidation = new Modification(35, new Composition.Composition(0, 0, 0, 1, 0), "Oxidation");
        public static readonly Modification DiMethylation = new Modification(36, new Composition.Composition(2, 4, 0, 0, 0), "Dimethyl");
        public static readonly Modification TriMethylation = new Modification(37, new Composition.Composition(3, 6, 0, 0, 0), "Trimethyl");

        public static readonly Modification Glutathione = new Modification(55, new Composition.Composition(10, 15, 3, 6, 1), "Glutathione");
        public static readonly Modification Cysteinyl = new Modification(312, new Composition.Composition(3, 5, 1, 2, 1), "Cysteinyl");
        public static readonly Modification Dehydro = new Modification(374, new Composition.Composition(0, -1, 0, 0, 0), "Dehydro");

        public static readonly Modification Itraq4Plex = new Modification(214, Data.Composition.Composition.Parse("H(12) C(4) 13C(3) N 15N O"), "iTRAQ4plex");
        public static readonly Modification Tmt6Plex = new Modification(737, Data.Composition.Composition.Parse("H(20) C(8) 13C(4) N 15N O(2)"), "TMT6plex");

        public static readonly Modification Nethylmaleimide = new Modification(108, Data.Composition.Composition.Parse("H(7) C(6) N O(2)"), "Nethylmaleimide");
        public static readonly Modification Nitrosyl = new Modification(275, Data.Composition.Composition.Parse("H(-1) N O"), "Nitrosyl");

        private static readonly Dictionary<string, Modification> NameToModMap;
        private static readonly Dictionary<string, IList<Modification>> MassToModMap;

        private static readonly Modification[] CommonModifications =
            {
                Acetylation,
                Carbamidomethylation,
                Carbamylation,
                Carboxymethylation,
                Deamidation,
                NipCam,
                Phosphorylation,
                PyroCarbamidomethyl,
                PyroGluE,
                PyroGluQ,
                Methylation,
                Oxidation,
                Itraq4Plex,
                Tmt6Plex,
                Glutathione,
                Dehydro,
                Nethylmaleimide,
                Nitrosyl
            };

        static Modification()
        {
            NameToModMap = new Dictionary<string, Modification>();
            MassToModMap = new Dictionary<string, IList<Modification>>();
            foreach (var modification in CommonModifications)
            {
                NameToModMap.Add(modification.Name, modification);
                var massStr = string.Format("{0:N3}", modification.Composition.Mass);
                //Console.WriteLine("{0} {1}", modification.Name, massStr);
                IList<Modification> modList;
                if (!MassToModMap.TryGetValue(massStr, out modList))
                {
                    modList = new List<Modification> {modification};
                    MassToModMap[massStr] = modList;
                }
                modList.Add(modification);
            }
        }
    }
}
