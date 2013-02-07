using System;
using System.Collections.Generic;
using System.IO;

namespace InformedProteomics.Backend.Scoring
{
    public class Misc
    {

        public static IEnumerable<string> GetPeptidesFromProtein(string protein, bool fullyTryptic,
                                                              int missedCleavageNumber)
        {
            var peptides = new HashSet<string>();

            for (var sp = 0; sp < protein.Length-4; sp++)
            {
                var numKR = 0;
                var lterm = sp == 0 || protein[sp - 1] == 'K' || protein[sp - 1] == 'R';

                for (var ep = sp; ep < Math.Min(sp + 50, protein.Length); ep++)
                {
                    var cterm = ep == protein.Length - 1 || protein[ep] == 'K' || protein[ep] == 'R';
                    if (cterm) numKR++;

                    if (ep - sp + 1 >=4 && numKR - 1 <= missedCleavageNumber &&
                        ((fullyTryptic && lterm && cterm) || !fullyTryptic && (lterm || cterm)))
                    {
                        peptides.Add(protein.Substring(sp, ep - sp + 1));
                    }
                }
            }

                return peptides;
        }

        public static HashSet<string> GetPeptidesFromFasta(string file, bool fullyTryptic, int missedCleavageNumber, bool reverse)// very slow.. for test only
        {
            var peptides = new HashSet<string>();
            var proteins = new List<string>();
            var sr = new StreamReader(file);

            var protein = "";
            while (true)
            {
                var s = sr.ReadLine();
                if (s == null || s.StartsWith(">"))
                {
                    if (protein.Length > 0)
                    {
                        if (reverse)
                        {
                            var revProtein = "";
                            for (var i = 0; i < protein.Length; i++)
                            {
                                revProtein += protein[protein.Length-i-1];
                            }
                            protein = revProtein;
                        }
                        proteins.Add(protein);
                        protein = "";
                    }
                    if (s == null) break;
                }
                else protein += s;
            }
            sr.Close();

            foreach (var p in proteins)
            {
                foreach (var pep in GetPeptidesFromProtein(p, fullyTryptic, missedCleavageNumber))
                {
                    peptides.Add(pep);
                }
            }

            return peptides;
        } 
    }
}
