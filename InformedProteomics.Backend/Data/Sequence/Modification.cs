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

        public static readonly Modification NoModification = new Modification(0, new Composition(0, 0, 0, 0, 0), "No modification");
        public static readonly Modification Carbamidomethylation = new Modification(4, new Composition(2, 3, 1, 1, 0), "Carbamidomethyl");
    }
}
