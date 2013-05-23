using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMolecule
    {
        public Sequence(IEnumerable<AminoAcid> aaArr)
        {
            Composition composition = Composition.Zero;
            SequenceString = "";

            foreach (AminoAcid aa in aaArr)
            {
                Add(aa);
                composition += aa.Composition;
                SequenceString += aa.Residue; // added by Kyowon jeong
            }

            Composition = composition + Composition.H2O; // +H2PO added by Kyowon Jeong
        }

		/*public Sequence (Composition composition, string sequence)
		{
			this.Composition = composition;
			this.SequenceString = sequence;
		}*/

        // fixed by Kyowon so that both constructors have amino acid list.
        public Sequence(Composition composition, string sequence, AminoAcidSet aminoAcidSet)
        {
            Composition = composition;
            SequenceString = sequence;
            foreach(var residue in SequenceString)
            {
                Add(aminoAcidSet.GetAminoAcid(residue));
            }
        }

        public Sequence(string sequence, AminoAcidSet aminoAcidSet)
        {
            Composition composition = Composition.Zero;
            SequenceString = sequence;
            foreach (var residue in SequenceString)
            {
                var aa = aminoAcidSet.GetAminoAcid(residue);
                Add(aa);
                composition += aa.Composition;
            }
            Composition = composition;
        }

        public Composition Composition { get; private set; }
    	public string SequenceString { get; set; }

        public double GetMass()
        {
            return Composition.GetMass();
        }

        public Composition GetComposition()
        {
            return Composition;
        }

        //added by kyowon jeong
        public double GetMass(int from, int to)
        {
            return GetComposition(from, to).GetMass();
        }

        //added by kyowon jeong
        public Composition GetComposition(int from, int to)
        {
            from = Math.Max(from, 0);
            to = Math.Min(to, Count);
            var composition = Composition.Zero;
            for (var i = from; i < to; i++)
                composition += this[i].Composition;
            return composition;
        }

        public IEnumerable<Composition> GetPrefixCompositions()
        {
            var compositions = new Composition[Count];
            var prefixComposition = Composition.Zero;
            var index = -1;
            foreach (var aa in this)
            {
                compositions[++index] = (prefixComposition += aa.Composition);
            }
            return compositions;
        }

        public IEnumerable<Composition> GetSuffixCompositions()
        {
            var compositions = new Composition[Count];
            var suffixComposition = Composition.Zero;
            for(var index = Count-1; index >= 0; --index)
            {
                compositions[index] = (suffixComposition += this[index].Composition);
            }
            return compositions;
        }

        public Ion GetPrecursorIon(int charge)
        {
            return new Ion(Composition, charge);
        }

        public Dictionary<string, Ion> GetProductIons(IEnumerable<IonType> ionTypes)
        {
            var ionTypeArr = ionTypes as IonType[] ?? ionTypes.ToArray();

            var productIonMap = new Dictionary<string, Ion>();

            // prefix
            var index = 0;
            foreach (var prefixComposition in GetPrefixCompositions())
            {
                ++index;
                foreach (var ionType in ionTypeArr.Where(ionType => ionType.IsPrefixIon))
                {
                    productIonMap.Add("(" + ionType.Name + ")" + index, ionType.GetIon(prefixComposition));
                }
            }

            // suffix
            index = 0;
            foreach (var suffixComposition in GetSuffixCompositions())
            {
                ++index;
                foreach (var ionType in ionTypeArr.Where(ionType => !ionType.IsPrefixIon))
                {
                    productIonMap.Add("(" + ionType.Name + ")" + index, ionType.GetIon(suffixComposition));
                }
            }

            return productIonMap;
        }
    }
}
