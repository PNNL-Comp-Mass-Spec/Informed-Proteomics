using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using InformedProteomics.Backend.Data.Sequence;

    /// <summary>
    /// Ion Type
    /// </summary>
    public class IonType
    {
        /// <summary>
        /// Name of ion
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Mass of ion
        /// </summary>
        public double Mass => OffsetComposition.Mass;

        /// <summary>
        /// Offset composition of the ion
        /// </summary>
        public Composition.Composition OffsetComposition { get; }

        /// <summary>
        /// Ion charge
        /// </summary>
        public int Charge { get; }

        /// <summary>
        /// If the ion is a prefix ion
        /// </summary>
        public bool IsPrefixIon { get; }

        /// <summary>
        /// BaseIonType of ion
        /// </summary>
        public BaseIonType BaseIonType { get; }

        /// <summary>
        /// Ion neutral loss
        /// </summary>
        public NeutralLoss NeutralLoss { get; }

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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="name"></param>
        /// <param name="offsetComposition"></param>
        /// <param name="charge"></param>
        /// <param name="isPrefixIon"></param>
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

        /// <summary>
        /// Get the m/z of <paramref name="cutMass"/> + (the offset mass)
        /// </summary>
        /// <param name="cutMass"></param>
        /// <returns></returns>
        public double GetMz(double cutMass)
        {
            return (cutMass + _offsetMass) / Charge + Constants.Proton;
        }

        /// <summary>
        /// Get the m/z of <paramref name="prefixComposition"/> + (the offset composition)
        /// </summary>
        /// <param name="prefixComposition"></param>
        /// <returns></returns>
        public double GetMz(Composition.Composition prefixComposition)
        {
            Debug.Assert(prefixComposition != null, "prefixComposition must not be null");
            return GetMz(prefixComposition.Mass);
        }

        /// <summary>
        /// Get the Ion with the cutComposition added
        /// </summary>
        /// <param name="cutComposition"></param>
        /// <returns></returns>
        public Ion GetIon(Composition.Composition cutComposition)
        {
            return new Ion(cutComposition + OffsetComposition, Charge);
        }

        /// <summary>
        /// Get possible ions for <paramref name="cutComposition"/> and <paramref name="terminalResidue"/>
        /// </summary>
        /// <param name="cutComposition"></param>
        /// <param name="terminalResidue"></param>
        /// <returns></returns>
        public IEnumerable<Ion> GetPossibleIons(Composition.Composition cutComposition, AminoAcid terminalResidue)
        {
            return BaseIonType.GetPossibleCompositions(terminalResidue)
                       .Select(offsetComposition => cutComposition + offsetComposition - NeutralLoss.Composition)
                       .Select(comp => new Ion(comp, Charge));
        }

        /// <summary>
        /// Get possible ions for <paramref name="sequence"/>
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public IEnumerable<Ion> GetPossibleIons(Sequence sequence)
        {
            var cutComposition = sequence.Composition;
            var terminalResidue = IsPrefixIon ? sequence[sequence.Count - 1] : sequence[0];
            return BaseIonType.GetPossibleCompositions(terminalResidue)
                       .Select(offsetComposition => cutComposition + offsetComposition - NeutralLoss.Composition)
                       .Select(comp => new Ion(comp, Charge));
        }

        /// <inheritdoc />
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
            string chargeStr;
            if (Charge == 1) chargeStr = "";
            else
            {
                var builder = new StringBuilder(Charge);
                for (var i = 0; i < Charge; i++) builder.Append("+");
                chargeStr = builder.ToString();
            }
            return BaseIonType.Symbol + index + chargeStr + NeutralLoss.Name;
        }

        /// <summary>
        /// Parse an ion from string <paramref name="s"/>
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
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

        /// <summary>
        /// Check for equality
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        protected bool Equals(IonType other)
        {
            return Charge == other.Charge && Equals(OffsetComposition, other.OffsetComposition) && IsPrefixIon == other.IsPrefixIon;
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (obj == null) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((IonType)obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            unchecked
            {
                //var hashCode = (this.Name != null ? this.Name.GetHashCode() : 0);
                var hashCode = Charge;
                hashCode = (hashCode * 397) ^ (OffsetComposition != null ? OffsetComposition.GetHashCode() : 0);
                hashCode = (hashCode * 397) ^ IsPrefixIon.GetHashCode();
                return hashCode;
            }
        }

        //}
    }
}
