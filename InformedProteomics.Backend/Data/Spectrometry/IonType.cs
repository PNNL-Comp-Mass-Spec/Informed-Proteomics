using System;
using System.Diagnostics;
using System.Text;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using System.Collections.Generic;
    using System.Linq;

    using InformedProteomics.Backend.Data.Sequence;

    public class IonType
    {
        public string Name { get; private set; }
        public Composition.Composition OffsetComposition { get; private set; }
        public int Charge { get; private set; }
        public bool IsPrefixIon { get; private set; }
        public BaseIonType BaseIonType { get; private set; }
        public NeutralLoss NeutralLoss { get; private set; }

        private readonly double _offsetMass;    // duplication but stored for performance

        internal IonType(
            string name, 
            Composition.Composition offsetComposition,
            int charge, 
            BaseIonType baseIonType,
            NeutralLoss neutralLoss
            )
        {
            Name = name;
            _offsetMass = offsetComposition.Mass;
            OffsetComposition = offsetComposition;
            Charge = charge;
            IsPrefixIon = baseIonType.IsPrefix;
            BaseIonType = baseIonType;
            NeutralLoss = neutralLoss;
        }

        [Obsolete("IonType object should be generated through IonTypeFactory")]
        public IonType(string name, Composition.Composition offsetComposition,
                       int charge, bool isPrefixIon)
        {
            Name = name;
            OffsetComposition = offsetComposition;
            Charge = charge;
            IsPrefixIon = isPrefixIon;
            BaseIonType = null;
            _offsetMass = offsetComposition.Mass;
        }

        public double GetMz(double cutMass)
        {
            return (cutMass + _offsetMass) / Charge + Constants.Proton;
        }

        public double GetMz(Composition.Composition prefixComposition)
        {
            Debug.Assert(prefixComposition != null, "prefixComposition must not be null");
            return GetMz(prefixComposition.Mass);
        }

        public Ion GetIon(Composition.Composition cutComposition)
        {
            return new Ion(cutComposition + OffsetComposition, Charge);
        }

        public IEnumerable<Ion> GetPossibleIons(Composition.Composition cutComposition, AminoAcid terminalResidue)
        {
            return this.BaseIonType.GetPossibleCompositions(terminalResidue)
                       .Select(offsetComposition => cutComposition + offsetComposition - this.NeutralLoss.Composition)
                       .Select(comp => new Ion(comp, Charge));
        }

        public IEnumerable<Ion> GetPossibleIons(Sequence sequence)
        {
            var cutComposition = sequence.Composition;
            var terminalResidue = IsPrefixIon ? sequence[sequence.Count - 1] : sequence[0];
            return this.BaseIonType.GetPossibleCompositions(terminalResidue)
                       .Select(offsetComposition => cutComposition + offsetComposition - this.NeutralLoss.Composition)
                       .Select(comp => new Ion(comp, Charge));
        } 

        public override string ToString()
        {
            return Name + "," + OffsetComposition + "," + _offsetMass +
                   "," + Charge + "," + IsPrefixIon;
        }

		/// <summary>
		/// Returns ion name with ion index (e.g. y++4-H2O => charge 2 y4 - H2O)
		/// </summary>
		/// <param name="index"></param>
		/// <returns></returns>
		public string GetName(int index)
		{
			String chargeStr;
			if (Charge == 1) chargeStr = "";
			else
			{
				var builder = new StringBuilder(Charge);
				for (var i = 0; i < Charge; i++) builder.Append("+");
				chargeStr = builder.ToString();
			}
			return BaseIonType.Symbol + index + chargeStr + NeutralLoss.Name;
		}

        [Obsolete("IonType object should be generated through IonTypeFactory")]
        public static IonType Parse(string s)
        {
            var t = s.Split(',');
            if (t.Length < 5) return null;
            var name = t[0];
            var composition = Composition.Composition.Parse(t[1]);
            var charge = int.Parse(t[3]);
            var isPrefixIon = bool.Parse(t[4]);
            return new IonType(name, composition, charge, isPrefixIon);
        }

        //public override bool Equals(object obj)
        //{
        //    var type = obj as IonType;
        //    if (type != null)
        //    {
        //        var other = type;
        //        return other.IsPrefixIon.Equals(IsPrefixIon) && other.OffsetComposition.Equals(OffsetComposition) &&
        //               other.Charge == Charge;
        //    }
        //    return false;
        //}

        //public override int GetHashCode()
        //{
        //    //return IsPrefixIon.GetHashCode() + OffsetComposition.GetHashCode() + Charge.GetHashCode();

        protected bool Equals(IonType other)
        {
            return /*string.Equals(this.Name, other.Name) &&*/ Charge == other.Charge && Equals(this.OffsetComposition, other.OffsetComposition) && this.IsPrefixIon == other.IsPrefixIon;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((IonType)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                //var hashCode = (this.Name != null ? this.Name.GetHashCode() : 0);
                var hashCode = Charge;
                hashCode = (hashCode * 397) ^ (this.OffsetComposition != null ? this.OffsetComposition.GetHashCode() : 0);
                hashCode = (hashCode * 397) ^ this.IsPrefixIon.GetHashCode();
                return hashCode;
            }
        }

        //}
    }
  
}
