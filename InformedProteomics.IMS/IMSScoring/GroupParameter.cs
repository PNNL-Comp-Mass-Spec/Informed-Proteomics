using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.IMS.IMSScoring
{
    public class GroupParameter
    {
        public int MassIndex { get; private set; }
        public int LocationIndex { get; private set; }
        public int FlankingResidueIndex { get; private set; }
        public int Charge { get; private set; }
        
        private const int MaxCharge = 3;
        private const int MinCharge = 2;
        private static Dictionary<Tuple<char, char>, int> flankingIndexDictionary = new Dictionary<Tuple<char, char>, int>(); 

        public GroupParameter(Composition cutComposition, char nTermAA, char cTermAA, Ion precursorIon)
        {
            MassIndex = GetMassIndex(precursorIon.Composition);
            Charge = Math.Max(MinCharge, Math.Min(precursorIon.Charge, MaxCharge));
            LocationIndex = GetLocationIndex(precursorIon.Composition, cutComposition);
            FlankingResidueIndex =  GetFlankingResidueIndex(nTermAA, cTermAA); 
        }

        public GroupParameter(Ion precursorIon) : this(precursorIon.Composition, ' ', ' ', precursorIon)
        {
            
        }

        internal GroupParameter(int massIndex, int locationIndex, int flankingResidueIndex, int charge)
        {
            MassIndex = massIndex;
            LocationIndex = locationIndex;
            FlankingResidueIndex = flankingResidueIndex;
            Charge = charge;
        }

        public GroupParameter GetPrecursorGroupParameter()
        {
            return new GroupParameter(MassIndex, int.MinValue, int.MinValue, Charge);
        } 

        public static GroupParameter Parse(string s)
        {
            var token = s.Split(' ');
            return token.Length != 4 ? null : new GroupParameter(int.Parse(token[0]), int.Parse(token[1]), int.Parse(token[2]), int.Parse(token[3]));
        }

        public static List<GroupParameter> GetAllFragmentGroupParameters(int maxCharge)
        {
            var ret = new List<GroupParameter>();
            for (var i = 1; i <= 3; i++)
            {
                for (var j = 1; j <= 4; j++)
                {
                    for (var k = 0; k <= 15; k++) // 16 for precursor
                    {
                        for (var c = 1; c <= maxCharge;c++ ){
                            ret.Add(new GroupParameter(i, j, k, c));
                        }
                    }
                }
            }
            return ret;
        }

        public static List<GroupParameter> GetAllPrecursorGroupParameters(int maxCharge)
        {
            var ret = new List<GroupParameter>();
            for (var i = 1; i <= 3; i++)
            {
                for (var j = 1; j <= 4; j++)
                {
                    for (var c = 1; c <= maxCharge; c++)
                    {
                        ret.Add(new GroupParameter(i, j, 16, c));
                    }
                }
            }
            return ret;
        }

        public override bool Equals(object o)
        {
            if (this == o) return true;
            var other = o as GroupParameter;
            if (other != null)
            {
                return
                    MassIndex == other.MassIndex && LocationIndex == other.LocationIndex && FlankingResidueIndex == other.FlankingResidueIndex && Charge == other.Charge;
            }
            return false;
        }

        public override int GetHashCode()
        {
            return MassIndex.GetHashCode() * LocationIndex.GetHashCode() * FlankingResidueIndex.GetHashCode() * Charge.GetHashCode();
        }

        public override string ToString()
        {
            return MassIndex + " " + LocationIndex + " " + FlankingResidueIndex + " " + Charge;
        }

        private static int GetMassIndex(Composition precursorComposition)
        {
            var m = precursorComposition.Mass;
            return m < 1200 ? 1 : (m < 2400 ? 2 : 3);
        }

        private static int GetLocationIndex(Composition precursorComposition, Composition cutComposition)
        {
            return (int)(cutComposition.Mass / precursorComposition.Mass * 4 + 1);
        }

        private static int GetResidueIndex(char residue)
        {
            const string residues = "ARNDCQEGHILKMFPSTWYV";
            return residues.IndexOf(residue);
        }

        private static int GetFlankingResidueIndex(char nTermAA, char cTermAA)
        {
            if (nTermAA.Equals(' ')) return 16;
            var key = new Tuple<char, char>(nTermAA, cTermAA);
            if (flankingIndexDictionary.ContainsKey(key))
                return flankingIndexDictionary[key];
            var index = ResidueTable[GetResidueIndex(nTermAA), GetResidueIndex(cTermAA)];
            flankingIndexDictionary[key] = index;
            return index;
        }

       

        static private readonly int[,] ResidueTable =
        {
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {14, 5, 15, 9, 14, 14, 9, 4, 7, 11, 11, 5, 14, 14, 2, 13, 13, 14, 14, 11},
            {8, 5, 8, 9, 8, 8, 8, 4, 7, 8, 8, 5, 8, 8, 2, 8, 8, 8, 8, 8},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {8, 5, 8, 9, 8, 8, 9, 4, 7, 8, 8, 5, 8, 8, 2, 8, 8, 8, 8, 8},
            {3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3},
            {6, 5, 6, 6, 6, 6, 6, 4, 7, 6, 6, 5, 6, 6, 2, 6, 6, 6, 6, 6},
            {10, 5, 10, 9, 10, 10, 9, 4, 7, 11, 10, 5, 10, 10, 2, 10, 10, 10, 10, 10},
            {10, 5, 10, 9, 10, 10, 9, 4, 7, 11, 11, 5, 10, 10, 2, 10, 10, 10, 10, 10},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1},
            {12, 5, 12, 9, 12, 12, 9, 4, 7, 11, 11, 5, 12, 12, 2, 13, 12, 12, 12, 11},
            {12, 5, 12, 9, 12, 12, 9, 4, 7, 11, 11, 5, 12, 12, 2, 13, 13, 12, 12, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {0, 5, 15, 9, 0, 0, 9, 4, 7, 11, 11, 5, 0, 0, 2, 13, 13, 0, 0, 11},
            {10, 5, 10, 9, 10, 10, 9, 4, 7, 11, 11, 5, 10, 10, 2, 10, 10, 10, 10, 11}
        };

    }
}
