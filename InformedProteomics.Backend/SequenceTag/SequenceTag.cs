using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTag : List<SequenceTagGraphEdge>, IEquatable<SequenceTag>, IComparable<SequenceTag>
    {
        public string HashString { get; private set; }
        public double Score;
        public double TotalMass { get; private set; }
        public double FrankingMass { get; private set; }
        public SequenceTag(IEnumerable<SequenceTagGraphEdge> edges)
            : base(edges)
        {
            var nominalMassArray = new byte[Count];
            for (var i = 0; i < Count; i++)
            {
                nominalMassArray[i] = this[i].NominalMass;
                TotalMass += this[i].Mass;
            }
            HashString = BitConverter.ToString(nominalMassArray);
            Score = 0;
        }

        public bool Equals(SequenceTag other)
        {
            return HashString.Equals(other.HashString);
        }

        public override int GetHashCode()
        {
            return HashString.GetHashCode();
        }

        public int CompareTo(SequenceTag other)
        {
            //return -Score.CompareTo(other.Score); // negative log
            return Score.CompareTo(other.Score);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((SequenceTag)obj);
        }

        public string[] GetTagStrings(out double[] rmse)
        {
            var listOfAminoAcids = new List<List<KeyValuePair<double, char>>>(Count);

            for (var j = 0; j < Count; j++)
            {
                var temp = this[j].AminoAcidList.Select(aa => new KeyValuePair<double, char>(aa.Value, aa.Key.Residue)).ToList();
                listOfAminoAcids.Add(temp);
            }

            var indexArray = new int[Count];
            var totalCombinations = listOfAminoAcids.Aggregate(1, (current, x) => current * x.Count);
            var result = new string[totalCombinations];
            rmse = new double[totalCombinations];

            for (var e = 0; e < totalCombinations; e ++)
            {
                var sb = new StringBuilder();
                for (var i = 0; i < indexArray.Length; i++)
                {
                    sb.Append(listOfAminoAcids[i][indexArray[i]].Value);
                    rmse[e] += listOfAminoAcids[i][indexArray[i]].Key * listOfAminoAcids[i][indexArray[i]].Key;
                }

                //increase indexArray
                for (var i = indexArray.Length - 1; i >= 0; i--)
                {
                    if (indexArray[i] == listOfAminoAcids[i].Count - 1) //reached the last item
                    {
                        indexArray[i] = 0;
                    }
                    else
                    {
                        indexArray[i]++;
                        break;
                    }
                }


                result[e] = sb.ToString();
                rmse[e] = Math.Sqrt(rmse[e] / Count);
            }

            Array.Sort(rmse, result);
            return result;
        }

        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

    }
}
