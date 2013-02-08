using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Scoring
{
    internal class FragmentSpectrumParameter : FragmentParameter
    {
        private readonly int _flankingResidueIndex;

        public FragmentSpectrumParameter(Sequence peptide, int residueNumber) : base(peptide, residueNumber)
        {
            _flankingResidueIndex = GetFlankingResidueIndex(peptide, residueNumber);  
        }

        private FragmentSpectrumParameter(int massIndex, int locationIndex, int flankingResidueIndex): base(massIndex, locationIndex)
        {
            _flankingResidueIndex = flankingResidueIndex;
        }

        private static int GetResidueIndex(char residue)
        {
            const string residues = "ARNDCQEGHILKMFPSTWYV";
            return residues.IndexOf(residue);
        }

        public new static FragmentSpectrumParameter Parse(string s)
        {
            var token = s.Split(' ');
            return token.Length != 3 ? null : new FragmentSpectrumParameter(int.Parse(token[0]), int.Parse(token[1]), int.Parse(token[2]));
        }

        public static List<FragmentSpectrumParameter> GetAllFragmentSpectrumParameters()
        {
            var ret = new List<FragmentSpectrumParameter>();
            for (var i = 1; i <= 3; i++)
            {
                for (var j = 1; j <= 4; j++)
                {
                    for (var k = 0; k <= 15; k++)
                    {
                        ret.Add(new FragmentSpectrumParameter(i, j, k));
                    }
                }
            }
            return ret;
        }

        private static int GetFlankingResidueIndex(Sequence peptide, int residueNumber)
        {
            var nterm = peptide[residueNumber - 1].Residue;
            var cterm = peptide[residueNumber].Residue;
            return ResidueTable[GetResidueIndex(nterm), GetResidueIndex(cterm)];
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

        public override bool Equals(object o)
        {
            if (this == o) return true;
            var other = o as FragmentSpectrumParameter;
            
            if (other != null)
            {
                return base.Equals(o) && _flankingResidueIndex == other._flankingResidueIndex;
            }
            return false;
        }

        public override int GetHashCode()
        {
            return base.GetHashCode() * _flankingResidueIndex.GetHashCode();
        }

        public new string ToFileString()
        {
            return base.ToFileString() + " " + _flankingResidueIndex;
        }
    }
}
