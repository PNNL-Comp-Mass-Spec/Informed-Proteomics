using System.Collections.Generic;
using System.IO;
using System.Xml.Linq;

namespace InformedProteomics.Backend.Data.Biology
{
    using System;

    [Serializable]
    public class Atom : IMatter, IEquatable<Atom>
    {
        public Atom(string code, double mass, int nominalMass, string name)
        {
            Code = code;
            Mass = mass;
            NominalMass = nominalMass;
            Name = name;
        }

        public Atom()
        {

        }

        public string Code { get; set; }

        public string Name { get; set; }

        public double Mass { get; set; }

        public int NominalMass { get; set; }

        public static readonly Atom[] AtomArr =
            {
                // Use Unimod values
                new Atom("H", H, 1, "Hydrogen"),
                new Atom("2H", 2.014101779, 2, "Deuterium"),
                new Atom("D", 2.014101779, 2, "Deuterium"),
                new Atom("Li", 7.016003, 7, "Lithium"),
                new Atom("C", C, 12, "Carbon"),
                new Atom("13C", Constants.C13, 13, "Carbon13"),
                new Atom("N", N, 14, "Nitrogen"),
                new Atom("15N", 15.00010897, 15, "Nitrogen15"),
                new Atom("O", O, 16, "Oxygen"),
                new Atom("18O", 17.9991603, 18, "Oxygen18"),
                new Atom("F", 18.99840322, 19, "Fluorine"),
                new Atom("Na", 22.9897677, 23, "Sodium"),
                new Atom("P", 30.973762, 31, "Phosphorous"),
                new Atom("S", S, 32, "Sulfur"),
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
                new Atom("Hg", 201.970617, 202, "Mercury"),
                new Atom("B", 11.00930554, 11, "Boron"),
                new Atom("As", 74.92159, 75, "Arsenic"),
                new Atom("Mg",  23.985043, 24, "Magnesium"),

                // Compounds
                // TODO: Move to somewhere more appropriate
                new Atom("Hex", 162.052824, 162, "Hexose"),
                new Atom("HexNAc", 203.079373, 203, "N-Acetylhexosamine"),
                new Atom("dHex", 146.057909, 146, "Fucose"),
                new Atom("NeuAc", 291.095417, 291, "N-acetyl neuraminic acid"),
                new Atom("NeuGc", 307.090331, 307, "N-glycoyl neuraminic acid"),
                new Atom("Hep", 192.063388, 192, "Heptose"),
                new Atom("Pent", 132.042257, 85, "Pentose"),
            };

        public const double C = 12.0;
        public const double H = 1.007825035;
        public const double N = 14.003074;
        public const double O = 15.99491463;
        public const double S = 31.9720707;

        public static Atom Get(string code)
        {
            return AtomMap[code];
        }

        private static readonly Dictionary<string, Atom> AtomMap;
        static Atom()
        {
            AtomMap = new Dictionary<string, Atom>();
            foreach (var atom in AtomArr)
            {
                AtomMap.Add(atom.Code,atom);
            }
        }

        private static void ReadFromPnnlOmicsElementDataXmlFile()
        {
            const string xmlFileName = @"..\..\..\PNNLOmicsElementData.xml";
            var xdocument = XDocument.Load(xmlFileName);

            var parameterBaseElement = xdocument.Element("parameters");
            if (parameterBaseElement == null)
            {
                throw new IOException("Problem reading xml file " + xmlFileName + "; Expected element 'parameters' but it was not found");
            }

            var elementIsotopes = parameterBaseElement.Element("ElementIsotopes");
            if (elementIsotopes == null)
            {
                throw new IOException("Problem reading xml file " + xmlFileName + "; Expected element 'ElementIsotopes' but it was not found");
            }

            var elements = elementIsotopes.Elements("Element");
            if (elements == null)
            {
                throw new IOException("Problem reading xml file " + xmlFileName + "; Expected element 'Element' but it was not found");
            }

            foreach (var element in elements)
            {
                var symbol = element.Element("Symbol").Value;
                var name = element.Element("Name").Value;
                var numIsotopes = int.Parse(element.Element("NumIsotopes").Value);
            }
        }

        public bool Equals(Atom other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return this.Name.Equals(other.Name);
        }

        public override bool Equals(object other) {
            return Equals(other as Atom);
        }

        public override int GetHashCode() {
            var result = 29;
            result = result * 13 + (Name == null ? 0 : Name.GetHashCode());
            result = result * 13 + (Code == null ? 0 : Code.GetHashCode());
            return result;
        }
    }
}
