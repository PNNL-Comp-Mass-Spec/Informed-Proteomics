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
            var composition = Data.Composition.Composition.Zero;

            foreach (var aa in aaArr)
            {
                Add(aa);
                composition += aa.Composition;
            }

            Composition = composition;
        }

        public Sequence(string sequence, AminoAcidSet aminoAcidSet)
        {
            var composition = Data.Composition.Composition.Zero;
            foreach (var residue in sequence)
            {
                var aa = aminoAcidSet.GetAminoAcid(residue);
                Add(aa);
                composition += aa.Composition;
            }
            Composition = composition;
        }

        public Composition.Composition Composition { get; private set; }

        public double GetMass()
        {
            return Composition.Mass;
        }

        public Composition.Composition GetComposition()
        {
            return Composition;
        }

        public double GetMass(int from, int to)
        {
            return GetComposition(from, to).Mass;
        }

        // from: inclusive
        // to: exclusive
        public Composition.Composition GetComposition(int from, int to)
        {
            from = Math.Max(from, 0);
            to = Math.Min(to, Count);
            var composition = Data.Composition.Composition.Zero;
            for (var i = from; i < to; i++)
                composition += this[i].Composition;
            return composition;
        }

        public IEnumerable<Composition.Composition> GetPrefixCompositions()
        {
            var compositions = new Composition.Composition[Count];
            var prefixComposition = Data.Composition.Composition.Zero;
            var index = -1;
            foreach (var aa in this)
            {
                compositions[++index] = (prefixComposition += aa.Composition);
            }
            return compositions;
        }

        public IEnumerable<Composition.Composition> GetSuffixCompositions()
        {
            var compositions = new Composition.Composition[Count];
            var suffixComposition = Data.Composition.Composition.Zero;
            for(var index = 0; index < Count; ++index)
            {
                compositions[index] = (suffixComposition += this[Count-1-index].Composition);
            }
            return compositions;
        }

        public Ion GetPrecursorIon(int charge)
        {
            return new Ion(Composition + Data.Composition.Composition.H2O, charge);
        }

        public Dictionary<Tuple<IonType,int>, Ion> GetProductIons(IEnumerable<IonType> ionTypes)
        {
            var ionTypeArr = ionTypes as IonType[] ?? ionTypes.ToArray();

            var productIonMap = new Dictionary<Tuple<IonType, int>, Ion>();

            // prefix
            foreach (var ionType in ionTypeArr.Where(ionType => ionType.IsPrefixIon))
            {
                var index = 0;
                foreach (var prefixComposition in GetPrefixCompositions())
                {
                    ++index;
                    productIonMap.Add(new Tuple<IonType,int> (ionType, index), ionType.GetIon(prefixComposition));
                }
            }

            // suffix
            foreach (var ionType in ionTypeArr.Where(ionType => !ionType.IsPrefixIon))
            {
                var index = 0;
                foreach (var suffixComposition in GetSuffixCompositions())
                {
                    ++index;
                    productIonMap.Add(new Tuple<IonType, int>(ionType, index), ionType.GetIon(suffixComposition));
                }
            }

            return productIonMap;
        }
    }
}
