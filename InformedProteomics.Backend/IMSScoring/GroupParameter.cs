using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.IMSScoring
{
    public class GroupParameter
    {
        public int MassIndex { get; private set; }
        public int LocationIndex { get; private set; }
        public int FlankingResidueIndex { get; private set; }
        public int Charge { get; private set; }

        public GroupParameter(Composition cutComposition, char nTermAA, char cTermAA, Ion precursorIon)
        {
            MassIndex = GetMassIndex(precursorIon.Composition);
            Charge = precursorIon.Charge;
            LocationIndex = GetLocationIndex(precursorIon.Composition, cutComposition);
            FlankingResidueIndex = GetFlankingResidueIndex(nTermAA, cTermAA); 
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

        public static List<GroupParameter> GetAllFragmentParameters(int maxCharge)
        {
            var ret = new List<GroupParameter>();
            for (var i = 1; i <= 3; i++)
            {
                for (var j = 1; j <= 4; j++)
                {
                    for (var k = 0; k <= 15; k++)
                    {
                        for (var c = 1; c <= maxCharge;c++ ){
                            ret.Add(new GroupParameter(i, j, k, c));
                        }
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

        public string ToFileString()
        {
            return MassIndex + " " + LocationIndex + " " + FlankingResidueIndex + " " + Charge;
        }

        private static int GetMassIndex(Composition precursorComposition)
        {
            var m = precursorComposition.GetMass();
            return m < 1200 ? 1 : (m < 2400 ? 2 : 3);
        }

        private static int GetLocationIndex(Composition precursorComposition, Composition cutComposition)
        {
            return (int)(cutComposition.GetMass() / precursorComposition.GetMass() * 4 + 1);
        }

        private static int GetResidueIndex(char residue)
        {
            const string residues = "ARNDCQEGHILKMFPSTWYV";
            return residues.IndexOf(residue);
        }

        private static int GetFlankingResidueIndex(char nTermAA, char cTermAA)
        {
            return ResidueTable[GetResidueIndex(nTermAA), GetResidueIndex(cTermAA)];
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
