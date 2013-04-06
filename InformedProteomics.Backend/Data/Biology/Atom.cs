using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Biology
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
                new Atom("H", 1.007825035, 1, "Hydrogen"),
                new Atom("2H", 2.014101779, 2, "Deuterium"),
                new Atom("Li", 7.016003, 7, "Lithium"),
                new Atom("C", 12.0, 12, "Carbon"),
                new Atom("13C", 13.00335483, 13, "Carbon13"),
                new Atom("N", 14.003074, 14, "Nitrogen"),
                new Atom("15N", 15.00010897, 15, "Nitrogen15"),
                new Atom("O", 15.99491463, 16, "Oxigen"),
                new Atom("18O", 17.9991603, 18, "Oxigen"),
                new Atom("F", 18.99840322, 19, "Fluorine"),
                new Atom("Na", 22.9897677, 23, "Sodium"),
                new Atom("P", 30.973762, 13, "Phosphorous"),
                new Atom("S", 31.9720707, 32, "Sulfur"),
                new Atom("Cl", 34.96885272, 35, "Chlorine"),
                new Atom("K", 38.9637074, 39, "Potassium"),
                new Atom("Ca", 39.9625906, 40, "Calcium"),
                new Atom("Fe", 55.9349393, 56, "Iron"),
                new Atom("Ni", 57.9353462, 58, "Nickel"),
                new Atom("Cu", 62.9295989, 63, "Copper"),
                new Atom("Zn", 63.9291448, 64, "Zinc"),
                new Atom("Br", 78.9183361, 79, "Bromine"),
                new Atom("Se", 79.9165196, 80, "Selenium"),
                new Atom("Mo", 97.9054073, 98, "Molybdenum"),
                new Atom("Ag", 106.905092, 107, "Silver"),
                new Atom("I", 126.904473, 127, "Iodine"),
                new Atom("Au", 196.966543, 197, "Gold"),
                new Atom("Hg", 201.970617, 202, "Mercury")
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
