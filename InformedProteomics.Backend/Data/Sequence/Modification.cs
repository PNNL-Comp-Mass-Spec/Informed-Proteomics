using InformedProteomics.Backend.Data.Science;
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

        public static readonly Modification Carbamidomethylation = new Modification(4, new Composition(2, 3, 1, 1, 0), "Carbamidomethyl");
    }
}
