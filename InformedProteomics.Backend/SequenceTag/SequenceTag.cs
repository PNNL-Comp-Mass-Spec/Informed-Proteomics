using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTag : List<SequenceTagGraphEdge>, IEquatable<SequenceTag>, IComparable<SequenceTag>
    {
        public string HashString { get; private set; }
        public double Score;
        public double TotalMass { get; private set; }
        public double FrankingMass { get; private set; }
        private readonly IList<DeconvolutedPeak> _peaks;
        private readonly Tolerance _tolerance;
        public HashSet<string> TagStrings;

        public SequenceTag(IEnumerable<SequenceTagGraphEdge> edges, IList<DeconvolutedPeak> peaks, Tolerance tolerance)
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
            _peaks = peaks;
            //TagStrings = null;
            _tolerance = tolerance;
            GetTagStrings();
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

        public void Merge(HashSet<string> tagStrings)
        {
            //if (TagStrings == null) return;
            foreach (var str in tagStrings) TagStrings.Add(str);
        }

        public HashSet<string> GetTagStrings()
        {
            if (TagStrings != null) return TagStrings;

            var listOfAminoAcids = new List<List<AminoAcid>>(Count);

            for (var j = 0; j < Count; j++)
            {
                listOfAminoAcids.Add(this[j].AminoAcidList);
            }

            var indexArray = new int[Count];
            var totalCombinations = listOfAminoAcids.Aggregate(1, (current, x) => current * x.Count);
            TagStrings = new HashSet<string>();
            //var result = new string[totalCombinations];
            //var rmse = new double[totalCombinations];
            var massTh = _tolerance.GetToleranceAsTh(_peaks[this[0].Node1].Mass);

            for (var e = 0; e < totalCombinations; e++)
            {
                var sb = new StringBuilder();
                var mass = 0d;
                for (var i = 0; i < indexArray.Length; i++)
                {
                    sb.Append(listOfAminoAcids[i][indexArray[i]].Residue);
                    mass += listOfAminoAcids[i][indexArray[i]].Mass;
                }

                var massGap = _peaks[this[indexArray.Length- 1].Node2].Mass - _peaks[this[0].Node1].Mass;
                var massError = Math.Abs(massGap - mass);

                if (massError < massTh) TagStrings.Add(sb.ToString());

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
                //result[e] = sb.ToString();
                //rmse[e] = Math.Sqrt(rmse[e] / Count);
            }
            return TagStrings;
        }
        /*
        public IList<string> GetTagStrings_old(Tolerance tolerance)
        {
            if (_tagString != null) return _tagString;
            
            //var listOfAminoAcids = new List<List<KeyValuePair<double, char>>>(Count);
            var listOfAminoAcids = new List<List<AminoAcid>>(Count);
            //var tolerance = new Tolerance(10);

            for (var j = 0; j < Count; j++)
            {
                //var temp = this[j].AminoAcidList.Select(aa => new KeyValuePair<double, char>(aa.Value, aa.Key.Residue)).ToList();
                listOfAminoAcids.Add(this[j].AminoAcidList);
            }

            var indexArray = new int[Count];
            var totalCombinations = listOfAminoAcids.Aggregate(1, (current, x) => current * x.Count);
            var result = new string[totalCombinations];
            var rmse = new double[totalCombinations];
            var massTh = tolerance.GetToleranceAsTh(_peaks[this[0].Node1].Mass);

            for (var e = 0; e < totalCombinations; e ++)
            {
                var sb = new StringBuilder();
                var mass = 0d;
                var bad = false;
                for (var i = 0; i < indexArray.Length; i++)
                {
                    sb.Append(listOfAminoAcids[i][indexArray[i]].Residue);
                    //rmse[e] += listOfAminoAcids[i][indexArray[i]].Key * listOfAminoAcids[i][indexArray[i]].Key;

                    mass += listOfAminoAcids[i][indexArray[i]].Mass;
                    var massGap = _peaks[this[i].Node2].Mass - _peaks[this[0].Node1].Mass;

                    var massError = Math.Abs(massGap - mass);
                    if (massError > massTh)
                    {
                        bad = true;
                        break;
                    }
                    rmse[e] = massError;
                }

                if (bad) rmse[e] = 10;

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
                //rmse[e] = Math.Sqrt(rmse[e] / Count);
            }

            Array.Sort(rmse, result);

            _tagString = new List<string>();
            for (var e = 0; e < totalCombinations; e++)
            {
                if (rmse[e] < 1) _tagString.Add(result[e]); 
            }
            return _tagString;
        }
        */
        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

    }
}
