using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagString : IEquatable<SequenceTagString>
    {
        public SequenceTagString(int scanNum, string sequence, bool isPrefix, double flankingMass, ActivationMethod actMethod = ActivationMethod.CID)
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
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;            
            return Equals(obj as SequenceTagString);
        }

        public bool Equals(SequenceTagString other)
        {
            if (other.ScanNum != ScanNum) return false;
            if (other.IsPrefix != IsPrefix) return false;
            if (!other.Sequence.Equals(Sequence)) return false;
            if (_mzComparer.GetBinNumber(FlankingMass) != _mzComparer.GetBinNumber(other.FlankingMass)) return false;
            return true;
        }

        private static readonly MzComparerWithBinning _mzComparer = new MzComparerWithBinning(28);

        public int ScanNum { get; private set; }
        public string Sequence { get; private set; }
        public bool IsPrefix { get; private set; }
        public double FlankingMass { get; private set; }

        public double TagMass
        {
            get
            {
                return _tagMass ??
                       (double)(_tagMass = AminoAcidSet.GetStandardAminoAcidSet().GetComposition(Sequence).Mass);
            }
        }

        public double? GetNTermFlankingMass(double? sequenceMass)
        {
            var baseIonTypes = _activationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;

            return GetNTermFlankingMass(sequenceMass, baseIonTypes[0], baseIonTypes[1]);
        }

        public double? GetNTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if (IsPrefix) return FlankingMass - prefixIonType.OffsetComposition.Mass;
            if (sequenceMass == null) return null;
            return sequenceMass - (FlankingMass - suffixIonType.OffsetComposition.Mass) - TagMass;
        }

        public double? GetCTermFlankingMass(double? sequenceMass)
        {
            var baseIonTypes = _activationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            return GetCTermFlankingMass(sequenceMass, baseIonTypes[0], baseIonTypes[1]);
        }

        public double? GetCTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if (!IsPrefix) return FlankingMass - suffixIonType.OffsetComposition.Mass;
            if (sequenceMass == null) return null;
            return sequenceMass - (FlankingMass - prefixIonType.OffsetComposition.Mass) - TagMass;
        }

        private double? _tagMass;
        private readonly ActivationMethod _activationMethod;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTagString()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

    }

}
