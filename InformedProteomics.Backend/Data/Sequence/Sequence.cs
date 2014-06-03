using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMolecule
    {
		public static readonly AminoAcidSet StandardAminoAcidSet;
		public static readonly AminoAcidSet StandardAminoAcidSetWithFixedCarbamidoMethyl;

	    static Sequence()
	    {
			StandardAminoAcidSet = new AminoAcidSet();
			StandardAminoAcidSetWithFixedCarbamidoMethyl = new AminoAcidSet(Modification.Carbamidomethylation);
	    }

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

        public IEnumerable<Cleavage> GetInternalCleavages()
        {
            var cleavages = new Cleavage[Count-1];
            var prefixComposition = Data.Composition.Composition.Zero;
            var suffixComposition = Data.Composition.Composition.Zero;
            for(var index = 0; index < Count-1; ++index)
            {
                cleavages[index] = new Cleavage(
                    prefixComposition += this[index].Composition,   // prefix
                    suffixComposition += this[Count - 1 - index].Composition    // suffix
                    );
            }
            return cleavages;
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

        // Parse peptide string from MS-GF+ results
        // e.g. +229.163C+57.021GLGGSGTPVDELDK+229.163C+57.021C+57.021QTHDNC+57.021YDQAK+229.163
        public static Sequence GetSequenceFromMsGfPlusPeptideStr(string msgfPlusPeptideStr)
        {
            const string aminoAcidRegex = @"[" + AminoAcid.StandardAminoAcidCharacters + "]";
            const string massRegex = @"[+-]?\d+\.\d+";

            if (!Regex.IsMatch(msgfPlusPeptideStr, "(" + aminoAcidRegex + "|" + massRegex + ")+")) return null;

            var stdAaSet = StandardAminoAcidSet;
            var aaList = new List<AminoAcid>();

            var matches = Regex.Matches(msgfPlusPeptideStr, "("+aminoAcidRegex+"|"+massRegex+")");
            AminoAcid aa = null;
            var mods = new List<Modification>();
            foreach (Match match in matches)
            {
                var element = match.Value;
                if (element.Length == 0) continue;
                if (element.Length == 1 && char.IsLetter(element[0]))   // amino acid
                {
                    if (aa != null)
                    {
                        aa = mods.Aggregate(aa, (current, mod) => new ModifiedAminoAcid(current, mod));
                        aaList.Add(aa);
                        mods = new List<Modification>();
                    }
                    aa = stdAaSet.GetAminoAcid(element[0]);
                    if(aa == null) throw new Exception("Unrecognized amino acid character: " + element[0]);
//                    Console.WriteLine("{0} {1} {2}", aa.Residue, aa.Composition, aa.GetMass());
                }
                else
                {
                    var modList = Modification.GetFromMass(element);
                    if (modList == null || modList.Count == 1) throw new Exception("Unrecognized modificaion mass: " + element);
                    var mod = modList[0];
                    mods.Add(mod);
//                    Console.WriteLine("{0} {1} {2}", mod.Name, mod.Composition, mod.Composition.AveragineMass);
                }
            }

            if (aa != null)
            {
                aa = mods.Aggregate(aa, (current, mod) => new ModifiedAminoAcid(current, mod));
                aaList.Add(aa);
            }

            return new Sequence(aaList);
        }

        //public static Sequence GetSequenceFromMsAlignSequenceStr(string msAlignSequenceStr)
        //{
        //    const string aminoAcidRegex = @"[" + AminoAcid.StandardAminoAcidCharacters + "]";
        //    const string massRegex = @"[+-]?\d+\.\d+";

        //    if (!Regex.IsMatch(msAlignSequenceStr, "(" + aminoAcidRegex + "|" + massRegex + ")+")) return null;

        //    var stdAaSet = StandardAminoAcidSet;
        //    var aaList = new List<AminoAcid>();

        //    var matches = Regex.Matches(msAlignSequenceStr, "(" + aminoAcidRegex + "|" + massRegex + ")");
        //    AminoAcid aa = null;
        //    var mods = new List<Modification>();
        //    foreach (Match match in matches)
        //    {
        //        var element = match.Value;
        //        if (element.Length == 0) continue;
        //        if (element.Length == 1 && char.IsLetter(element[0]))   // amino acid
        //        {
        //            if (aa != null)
        //            {
        //                aa = mods.Aggregate(aa, (current, mod) => new ModifiedAminoAcid(current, mod));
        //                aaList.Add(aa);
        //                mods = new List<Modification>();
        //            }
        //            aa = stdAaSet.GetAminoAcid(element[0]);
        //            if (aa == null) throw new Exception("Unrecognized amino acid character: " + element[0]);
        //            //                    Console.WriteLine("{0} {1} {2}", aa.Residue, aa.Composition, aa.GetMass());
        //        }
        //        else
        //        {
        //            var modList = Modification.GetFromMass(element);
        //            if (modList == null || modList.Count == 1) throw new Exception("Unrecognized modificaion mass: " + element);
        //            var mod = modList[0];
        //            mods.Add(mod);
        //            //                    Console.WriteLine("{0} {1} {2}", mod.Name, mod.Composition, mod.Composition.AveragineMass);
        //        }
        //    }

        //    if (aa != null)
        //    {
        //        aa = mods.Aggregate(aa, (current, mod) => new ModifiedAminoAcid(current, mod));
        //        aaList.Add(aa);
        //    }

        //    return new Sequence(aaList);
        //}
    }
}
