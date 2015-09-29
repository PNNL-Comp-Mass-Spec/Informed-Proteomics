using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMolecule
    {
        public Sequence(IEnumerable<AminoAcid> aaArr)
        {
            var aminoAcids = aaArr as IList<AminoAcid> ?? aaArr.ToList();

            PrefixComposition = new Composition.Composition[aminoAcids.Count + 1];
            PrefixComposition[0] = Data.Composition.Composition.Zero;
            for (var i = 0; i < aminoAcids.Count; i++)
            {
                var aa = aminoAcids[i] ?? AminoAcid.Empty;
                Add(aa);
                PrefixComposition[i+1] = PrefixComposition[i] + aa.Composition;
            }

            Composition = PrefixComposition[aminoAcids.Count];
            _prefixMass = PrefixComposition.Select(c => c.Mass).ToArray();
        }

        public Sequence(string sequence, AminoAcidSet aminoAcidSet): this(sequence.Select(aminoAcidSet.GetAminoAcid))
        {
        }

        public string GetModificationString()
        {
            var sb = new StringBuilder();
            for (var i = 0; i < Count; i++)
            {
                var modAa = this[i] as ModifiedAminoAcid;
                if (modAa == null) continue;

                if (sb.Length > 0) sb.Append(",");
                sb.AppendFormat("{0} {1}", modAa.Modification.Name, i);
            }
            return sb.ToString();
        }

        // modStr: E.g. Acetyl 0,Oxidation 1,Oxidation 20,Oxidation 27
        public static Sequence CreateSequence(string sequence, string modStr, AminoAcidSet aminoAcidSet)
        {
            if(string.IsNullOrEmpty(modStr)) return new Sequence(sequence, aminoAcidSet);

            var indexModMap = new Dictionary<int, Modification>();
            foreach (var modIns in modStr.Split(','))
            {
                var token = modIns.Split(' ');
                if (token.Length != 2) return null; // invalid modStr
                var mod = Modification.Get(token[0]);
                if (mod == null) return null;
                var index = Convert.ToInt32(token[1])-1;
                indexModMap.Add(index, mod);
            }

            var aaList = new List<AminoAcid>();
            
            for (var i=0; i<sequence.Length; i++)
            {
                var residue = sequence[i];
                var aa = aminoAcidSet.GetAminoAcid(residue);
                if (i == 0 && indexModMap.ContainsKey(-1))  // N-term modification
                {
                    var nTermMod = indexModMap[-1];
                    aa = new ModifiedAminoAcid(aa, nTermMod);
                }
                Modification mod;
                if (indexModMap.TryGetValue(i, out mod))
                {
                    var modifiedAa = new ModifiedAminoAcid(aa, mod);
                    aaList.Add(modifiedAa);
                }
                else
                {
                    aaList.Add(aa);
                }
            }

            return new Sequence(aaList);
        }

        // 1-based: PrefixComposition[0] = Composition.Zero
        public Composition.Composition[] PrefixComposition { get; private set; }
        public Composition.Composition Composition { get; private set; }
        public double Mass { get { return Composition.Mass; } }

        private readonly double[] _prefixMass;

        // from: inclusive
        // to: exclusive
        public double GetMass(int from, int to)
        {
            return _prefixMass[to] - _prefixMass[from];
        }

        // from: inclusive
        // to: exclusive
        public Composition.Composition GetComposition(int from, int to)
        {
            return PrefixComposition[to] - PrefixComposition[from];
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
            //var compositions = new Composition.Composition[Count];
            //var prefixComposition = Data.Composition.Composition.Zero;
            //var index = -1;
            //foreach (var aa in this)
            //{
            //    compositions[++index] = (prefixComposition += aa.Composition);
            //}
            //return compositions;
            return PrefixComposition.Skip(1);
        }

        public IEnumerable<Composition.Composition> GetSuffixCompositions()
        {
            //var compositions = new Composition.Composition[Count];
            //var suffixComposition = Data.Composition.Composition.Zero;
            //for(var index = 0; index < Count; ++index)
            //{
            //    compositions[index] = (suffixComposition += this[Count-1-index].Composition);
            //}
            //return compositions;
            return PrefixComposition.Reverse().Select(c => Composition - c).Take(Count);
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

            var stdAaSet = AminoAcidSet.GetStandardAminoAcidSet();
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
