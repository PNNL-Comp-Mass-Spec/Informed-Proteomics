using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Scoring
{
    internal class FragmentParameter
    {
        private readonly int _massIndex, _locationIndex;

        public FragmentParameter(Sequence peptide, int residueNumber)
        {
            _massIndex = GetMassIndex(peptide);
            _locationIndex = GetLocationIndex(peptide, residueNumber);
        }

        internal FragmentParameter(int massIndex, int locationIndex)
        {
            _massIndex = massIndex;
            _locationIndex = locationIndex;
        }

        public static FragmentParameter Parse(string s)
        {
            var token = s.Split(' ');
            return token.Length != 2 ? null : new FragmentParameter(int.Parse(token[0]), int.Parse(token[1]));
        }

        public static List<FragmentParameter> GetAllFragmentParameters()
        {
            var ret = new List<FragmentParameter>();
            for (var i = 1; i <= 3; i++)
            {
                for (var j = 1; j <= 4; j++)
                {
                    ret.Add(new FragmentParameter(i, j));
                }
            }
            return ret;
        }

        private static int GetMassIndex(Sequence peptide)
        {
            var len = peptide.Count;
            return len < 10 ? 1 : (len < 20 ? 2 : 3);
        }

        private static int GetLocationIndex(Sequence peptide, int residueNumber)
        {
            return (int) (peptide.GetMass(0, residueNumber)/peptide.GetMass()*4 + 1);
        }

        public override bool Equals(object o)
        {
            if (this == o) return true;
            var other = o as FragmentParameter;
            if (other != null)
            {
                return
                    _massIndex == other._massIndex && _locationIndex == other._locationIndex;
            }
            return false;
        }

        public override int GetHashCode()
        {
            return _massIndex.GetHashCode()*_locationIndex.GetHashCode();
        }

        public string ToFileString()
        {
            return _massIndex + " " + _locationIndex;
        }
    }
}
