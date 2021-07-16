using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// A sequence of Amino Acids, with modifications
    /// </summary>
    public class Sequence : List<AminoAcid>, IMolecule
    {
        // Ignore Spelling: acetyl

        /// <summary>
        /// Build a sequence from the supplied list of amino acids
        /// </summary>
        /// <param name="aaArr"></param>
        public Sequence(IEnumerable<AminoAcid> aaArr)
        {
            var aminoAcids = aaArr as IList<AminoAcid> ?? aaArr.ToList();

            PrefixComposition = new Composition.Composition[aminoAcids.Count + 1];
            PrefixComposition[0] = Data.Composition.Composition.Zero;
            for (var i = 0; i < aminoAcids.Count; i++)
            {
                var aa = aminoAcids[i] ?? AminoAcid.Empty;
                Add(aa);
                PrefixComposition[i + 1] = PrefixComposition[i] + aa.Composition;
            }

            Composition = PrefixComposition[aminoAcids.Count];
            _prefixMass = PrefixComposition.Select(c => c.Mass).ToArray();
        }

        /// <summary>
        /// Build a sequence from the supplied character sequence, using the provided amino acid set
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="aminoAcidSet"></param>
        public Sequence(string sequence, AminoAcidSet aminoAcidSet) : this(sequence.Select(aminoAcidSet.GetAminoAcid))
        {
        }

        /// <summary>
        /// Get the list (as a string) of the modifications and their locations in this sequence
        /// </summary>
        /// <returns>Human readable list of modifications</returns>
        public string GetModificationString()
        {
            var sb = new StringBuilder();
            for (var i = 0; i < Count; i++)
            {
                if (this[i] is not ModifiedAminoAcid modAa)
                {
                    continue;
                }

                if (sb.Length > 0)
                {
                    sb.Append(",");
                }

                sb.AppendFormat("{0} {1}", modAa.Modification.Name, i);
            }
            return sb.ToString();
        }

        /// <summary>
        /// Create a sequence using the supplied character sequence, modifications, and amino acid set
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="modStr">E.g. Acetyl 0,Oxidation 1,Oxidation 20,Oxidation 27</param>
        /// <param name="aminoAcidSet"></param>
        /// <returns>Sequence object</returns>
        public static Sequence CreateSequence(string sequence, string modStr, AminoAcidSet aminoAcidSet)
        {
            if (string.IsNullOrEmpty(modStr))
            {
                return new Sequence(sequence, aminoAcidSet);
            }

            var indexModMap = new Dictionary<int, Modification>();
            foreach (var modIns in modStr.Split(','))
            {
                var token = modIns.Split(' ');
                if (token.Length != 2)
                {
                    return null; // invalid modStr
                }

                var mod = Modification.Get(token[0]);
                if (mod == null)
                {
                    return null;
                }

                var index = Convert.ToInt32(token[1]) - 1;
                indexModMap.Add(index, mod);
            }

            var aaList = new List<AminoAcid>();

            for (var i = 0; i < sequence.Length; i++)
            {
                var residue = sequence[i];
                var aa = aminoAcidSet.GetAminoAcid(residue);
                if (i == 0 && indexModMap.ContainsKey(-1))  // N-term modification
                {
                    var nTermMod = indexModMap[-1];
                    aa = new ModifiedAminoAcid(aa, nTermMod);
                }

                if (indexModMap.TryGetValue(i, out var mod))
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

        /// <summary>
        /// 1-based: PrefixComposition[0] = Composition.Zero
        /// </summary>
        public Composition.Composition[] PrefixComposition { get; }

        /// <summary>
        /// Composition of the entire sequence
        /// </summary>
        public Composition.Composition Composition { get; }

        /// <summary>
        /// Mass of the entire sequence
        /// </summary>
        public double Mass => Composition.Mass;

        private readonly double[] _prefixMass;

        /// <summary>
        /// Get the mass from residue <paramref name="from"/> (inclusive) to residue <paramref name="to"/> (exclusive)
        /// </summary>
        /// <param name="from">inclusive</param>
        /// <param name="to">exclusive</param>
        /// <returns>Summed mass of the matching residues</returns>
        public double GetMass(int from, int to)
        {
            return _prefixMass[to] - _prefixMass[from];
        }

        /// <summary>
        /// Get the composition from residue <paramref name="from"/> (inclusive) to residue <paramref name="to"/> (exclusive)
        /// </summary>
        /// <param name="from">inclusive</param>
        /// <param name="to">exclusive</param>
        /// <returns>Composition of the matching residues</returns>
        public Composition.Composition GetComposition(int from, int to)
        {
            return PrefixComposition[to] - PrefixComposition[from];
        }

        /// <summary>
        /// Get the internal cleavages
        /// </summary>
        /// <returns>List of internal cleavage points</returns>
        public IEnumerable<Cleavage> GetInternalCleavages()
        {
            var cleavages = new Cleavage[Count - 1];
            var prefixComposition = Data.Composition.Composition.Zero;
            var suffixComposition = Data.Composition.Composition.Zero;
            for (var index = 0; index < Count - 1; ++index)
            {
                cleavages[index] = new Cleavage(
                    prefixComposition += this[index].Composition,   // prefix
                    this[index],
                    suffixComposition += this[Count - 1 - index].Composition,    // suffix
                    this[index + 1]
                    );
            }
            return cleavages;
        }

        /// <summary>
        /// Get the prefix compositions
        /// </summary>
        /// <returns>List of prefix compositions</returns>
        public IEnumerable<Composition.Composition> GetPrefixCompositions()
        {
            return PrefixComposition.Skip(1);
        }

        /// <summary>
        /// Get the suffix compositions
        /// </summary>
        /// <returns>List of suffix compositions</returns>
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

        /// <summary>
        /// Get the precursor ion
        /// </summary>
        /// <param name="charge"></param>
        /// <returns>Precursor ion for the given charge state</returns>
        public Ion GetPrecursorIon(int charge)
        {
            return new Ion(Composition + Data.Composition.Composition.H2O, charge);
        }

        /// <summary>
        /// Get the product ions of the specified types
        /// </summary>
        /// <param name="ionTypes"></param>
        /// <returns>Dictionary of product ions, where keys are Tuples of IonType and int, while are Ion objects</returns>
        public Dictionary<Tuple<IonType, int>, Ion> GetProductIons(IEnumerable<IonType> ionTypes)
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
                    productIonMap.Add(new Tuple<IonType, int>(ionType, index), ionType.GetIon(prefixComposition));
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

        /// <summary>
        /// Parse peptide string from MS-GF+ results
        /// </summary>
        /// <param name="msgfPlusPeptideStr">string of format "+229.163C+57.021GLGGSGTPVDELDK+229.163C+57.021C+57.021QTHDNC+57.021YDQAK+229.163"</param>
        /// <returns>Sequence object</returns>
        public static Sequence GetSequenceFromMsGfPlusPeptideStr(string msgfPlusPeptideStr)
        {
            const string aminoAcidRegex = "[" + AminoAcid.StandardAminoAcidCharacters + "]";
            const string massRegex = @"[+-]?\d+\.\d+";

            if (!Regex.IsMatch(msgfPlusPeptideStr, "(" + aminoAcidRegex + "|" + massRegex + ")+"))
            {
                return null;
            }

            var stdAaSet = AminoAcidSet.GetStandardAminoAcidSet();
            var aaList = new List<AminoAcid>();

            var matches = Regex.Matches(msgfPlusPeptideStr, "(" + aminoAcidRegex + "|" + massRegex + ")");
            AminoAcid aa = null;
            var mods = new List<Modification>();
            foreach (Match match in matches)
            {
                var element = match.Value;
                if (element.Length == 0)
                {
                    continue;
                }

                if (element.Length == 1 && char.IsLetter(element[0]))   // amino acid
                {
                    if (aa != null)
                    {
                        aa = mods.Aggregate(aa, (current, mod) => new ModifiedAminoAcid(current, mod));
                        aaList.Add(aa);
                        mods = new List<Modification>();
                    }
                    aa = stdAaSet.GetAminoAcid(element[0]);
                    if (aa == null)
                    {
                        throw new Exception("Unrecognized amino acid character: " + element[0]);
                    }
                    //                    Console.WriteLine("{0} {1} {2}", aa.Residue, aa.Composition, aa.GetMass());
                }
                else
                {
                    var modList = Modification.GetFromMass(element);
                    if (modList == null || modList.Count == 1)
                    {
                        throw new Exception("Unrecognized modification mass: " + element);
                    }

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
        //            if (modList == null || modList.Count == 1) throw new Exception("Unrecognized modification mass: " + element);
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
