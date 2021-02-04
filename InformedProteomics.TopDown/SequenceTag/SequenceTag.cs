using System;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.SequenceTag
{
    public class SequenceTag : IEquatable<SequenceTag>
    {
        public SequenceTag(int scanNum, string sequence, bool isPrefix, double flankingMass, ActivationMethod actMethod = ActivationMethod.CID)
        {
            ScanNum = scanNum;
            Sequence = sequence;
            IsPrefix = isPrefix;
            FlankingMass = flankingMass;
            _activationMethod = actMethod;
        }

        public override int GetHashCode()
        {
            var massBinIndex = _mzComparer.GetBinNumber(FlankingMass);
            var hashStr = string.Format("{0}{1}{2}", IsPrefix ? 1 : 0, massBinIndex, Sequence);
            return hashStr.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
            {
                return false;
            }

            if (ReferenceEquals(this, obj))
            {
                return true;
            }

            if (obj.GetType() != this.GetType())
            {
                return false;
            }

            return Equals(obj as SequenceTag);
        }

        public bool Equals(SequenceTag other)
        {
            if (other == null)
            {
                return false;
            }

            if (other.ScanNum != ScanNum)
            {
                return false;
            }

            if (other.IsPrefix != IsPrefix)
            {
                return false;
            }

            if (!other.Sequence.Equals(Sequence))
            {
                return false;
            }

            if (_mzComparer.GetBinNumber(FlankingMass) != _mzComparer.GetBinNumber(other.FlankingMass))
            {
                return false;
            }

            return true;
        }

        private static readonly MzComparerWithBinning _mzComparer = new MzComparerWithBinning(28);

        public int ScanNum { get; }
        public string Sequence { get; }
        public bool IsPrefix { get; }
        public double FlankingMass { get; }

        public double TagMass => _tagMass ?? (double)(_tagMass = AminoAcidSet.GetStandardAminoAcidSet().GetComposition(Sequence).Mass);

        public double? GetNTermFlankingMass(double? sequenceMass)
        {
            var baseIonTypes = _activationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;

            return GetNTermFlankingMass(sequenceMass, baseIonTypes[0], baseIonTypes[1]);
        }

        public double? GetNTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if (IsPrefix)
            {
                return FlankingMass - prefixIonType.OffsetComposition.Mass;
            }

            if (sequenceMass == null)
            {
                return null;
            }

            return sequenceMass - (FlankingMass - suffixIonType.OffsetComposition.Mass) - TagMass;
        }

        public double? GetCTermFlankingMass(double? sequenceMass)
        {
            var baseIonTypes = _activationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            return GetCTermFlankingMass(sequenceMass, baseIonTypes[0], baseIonTypes[1]);
        }

        public double? GetCTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if (!IsPrefix)
            {
                return FlankingMass - suffixIonType.OffsetComposition.Mass;
            }

            if (sequenceMass == null)
            {
                return null;
            }

            return sequenceMass - (FlankingMass - prefixIonType.OffsetComposition.Mass) - TagMass;
        }

        private double? _tagMass;
        private readonly ActivationMethod _activationMethod;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTag()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }
    }

    /*
    public class SequenceTag : List<SequenceTagGraphEdge>
    {
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
            //HashString = BitConverter.ToString(nominalMassArray);
            Score = 0;
            _peaks = peaks;
            //TagStrings = null;
            _tolerance = tolerance;
            GetTagStrings();
        }

        public void Merge(HashSet<string> tagStrings)
        {
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
            //var rootMeanSquareError = new double[totalCombinations];
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
                //rootMeanSquareError[e] = Math.Sqrt(rootMeanSquareError[e] / Count);
            }
            return TagStrings;
        }

        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
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
              if (obj == null) return false;
              if (ReferenceEquals(this, obj)) return true;
              if (obj.GetType() != this.GetType()) return false;
              return Equals((SequenceTag)obj);
          }
    }
    */
}
