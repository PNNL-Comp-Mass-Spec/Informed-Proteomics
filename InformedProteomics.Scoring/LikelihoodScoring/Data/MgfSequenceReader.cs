using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class MgfSequenceReader: ISequenceReader
    {
        private static readonly Dictionary<string, Tuple<AminoAcid, List<Modification>>> Modifications;

        private static readonly AminoAcidSet StandardAminoAcidSet;

        static MgfSequenceReader()
        {
            StandardAminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);
            Modifications = new Dictionary<string, Tuple<AminoAcid, List<Modification>>>();

            AddModification("99.032", 'G', Modification.Acetylation);
            AddModification("113.048", 'A', Modification.Acetylation);
            AddModification("129.043", 'S', Modification.Acetylation);
            AddModification("141.079", 'V', Modification.Acetylation);
            AddModification("143.059", 'T', Modification.Acetylation);
            AddModification("147.035", 'M', Modification.Oxidation);
            AddModification("157.038", 'D', Modification.Acetylation);
            AddModification("160.03", 'C', Modification.Carbamidomethylation);
            AddModification("171.054", 'E', Modification.Acetylation);
            AddModification("173.051", 'M', Modification.Acetylation);
            AddModification("189.046", 'F', Modification.Acetylation);
            AddModification("202.041", 'C', new List<Modification> { Modification.Carbamidomethylation, Modification.Acetylation });
        }

        private static void AddModification(string modMass, char residue, Modification modType)
        {
            var aminoAcid = StandardAminoAcidSet.GetAminoAcid(residue);
            var modList = new Tuple<AminoAcid, List<Modification>>(aminoAcid, new List<Modification> {modType});
            Modifications.Add(modMass, modList);
        }

        private static void AddModification(string modMass, char residue, List<Modification> modifications)
        {
            var aminoAcid = StandardAminoAcidSet.GetAminoAcid(residue);
            var modList = new Tuple<AminoAcid, List<Modification>>(aminoAcid, modifications);
            Modifications.Add(modMass, modList);
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
                        throw new Exception("Unrecognized modification mass: " + element);
                    }

//                    if (modList == null || modList.Count == 1) throw new Exception("Unrecognized modification mass: " + element);
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
