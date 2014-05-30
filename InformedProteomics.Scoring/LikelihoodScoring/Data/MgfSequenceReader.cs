using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class MgfSequenceReader: ISequenceReader
    {
        private readonly static Dictionary<string, Tuple<AminoAcid, List<Modification>>> Modifications;
        private readonly static AminoAcidSet StandardAminoAcidSet;
        static MgfSequenceReader()
        {
            StandardAminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);
            Modifications = new Dictionary<string, Tuple<AminoAcid, List<Modification>>>();
            Modifications.Add("99.032",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('G'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("113.048",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('A'), 
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("129.043",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('S'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("141.079",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('V'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("143.059",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('T'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("147.035",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('M'),
                                                        new List<Modification> { Modification.Oxidation }));
            Modifications.Add("157.038",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('D'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("160.03",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('C'),
                                                        new List<Modification> { Modification.Carbamidomethylation }));
            Modifications.Add("171.054",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('E'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("173.051",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('M'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("189.046",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('F'),
                                                        new List<Modification> { Modification.Acetylation }));
            Modifications.Add("202.041",
                new Tuple<AminoAcid, List<Modification>>(StandardAminoAcidSet.GetAminoAcid('C'),
                                                        new List<Modification> { Modification.Carbamidomethylation,
                                                                                 Modification.Acetylation }));
        }

        public Sequence GetSequence(string sequence)
        {
            const string aminoAcidRegex = @"[" + AminoAcid.StandardAminoAcidCharacters + "]";
            const string massRegex = @"\(\d+\.\d+\(";
            char[] parens = {'(', ')'};

            if (!Regex.IsMatch(sequence, "(" + aminoAcidRegex + "|" + massRegex + ")+")) return null;

            var stdAaSet = StandardAminoAcidSet;
            var aaList = new List<AminoAcid>();

            var matches = Regex.Matches(sequence, "(" + aminoAcidRegex + "|" + massRegex + ")");
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
                        mods.Clear();
                    }
                    aa = stdAaSet.GetAminoAcid(element[0]);
                    if (aa == null) throw new Exception("Unrecognized amino acid character: " + element[0]);
                    //                    Console.WriteLine("{0} {1} {2}", aa.Residue, aa.Composition, aa.GetMass());
                }
                else
                {
                    element = element.Trim(parens);
                    IList<Modification> modList;
                    AminoAcid modAa;
                    try
                    {
                        modList = Modifications[element].Item2;
                        modAa = Modifications[element].Item1;
                    }
                    catch (KeyNotFoundException)
                    {
                        throw new Exception("Unrecognized modificaion mass: " + element);
                    }

//                    if (modList == null || modList.Count == 1) throw new Exception("Unrecognized modificaion mass: " + element);
                    aa = modAa;
                    mods.AddRange(modList);
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
    }
}
