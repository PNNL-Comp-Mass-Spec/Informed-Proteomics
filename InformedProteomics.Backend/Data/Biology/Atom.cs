using System.Collections.Generic;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Science
{
    public class Atom : IMatter
    {
        public Atom(string code, double mass, int nominalMass, string name)
        {
            Code = code;
            Mass = mass;
            NominalMass = nominalMass;
            Name = name;
        }

        public string Name { get; private set; }

        public double Mass { get; private set; }

        public int NominalMass { get; private set; }

        public string Code { get; private set; }

        public double GetMass()
        {
            throw new System.NotImplementedException();
        }

        public static readonly Atom[] AtomArr =
            {
                new Atom("C", 12.0, 12, "Carbon"),
                new Atom("H", 1.007825035, 1, "Hydrogen"),
                new Atom("N", 14.003074, 14, "Nitrogen"),
                new Atom("O", 15.99491463, 16, "Oxigen"),
                new Atom("S", 31.9720707, 32, "Sulfur"),
            };

        public static Atom Get(string code)
        {
            return AtomMap[code];
        }

        private static readonly Dictionary<string,Atom> AtomMap = new Dictionary<string, Atom>();
        static Atom()
        {
            foreach (var atom in AtomArr)
            {
                AtomMap.Add(atom.Code,atom);
            }
        }
    }
}
