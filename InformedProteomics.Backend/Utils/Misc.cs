using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace InformedProteomics.Backend.Utils
{
    public class Misc
    {
        public static DateTime GetBuildDateFromVersion()
        {
            var version = Assembly.GetEntryAssembly().GetName().Version;
            return GetBuildDateFromVersion(version);
        }

        public static DateTime GetBuildDateFromVersion(Version version)
        {
            var buildDateTime = new DateTime(2000, 1, 1).Add(
                new TimeSpan(TimeSpan.TicksPerDay * version.Build + // days since 1 January 2000
                             TimeSpan.TicksPerSecond * 2 * version.Revision)); // seconds since midnight, (multiply by 2 to get original)

            return buildDateTime;
        }

        public static string GetBuildDateTextFromVersion()
        {
            var version = Assembly.GetEntryAssembly().GetName().Version;
            return GetBuildDateTextFromVersion(version);
        }

        public static string GetBuildDateTextFromVersion(Version version)
        {
            var buildDateTime = GetBuildDateFromVersion();
            return buildDateTime.ToString("MMMM d, yyyy");
        }

        public static List<string> GetPeptidesFromTxt(string fileName)
        {
            var peptides = new List<string>();

            var stremaReader = new StreamReader(fileName);
            string s;

            while ((s = stremaReader.ReadLine()) != null)
            {
                s = s.Replace("C!", "C");
                s = s.Substring(s.IndexOf('.') + 1, s.LastIndexOf('.') - s.IndexOf('.') - 1);
                peptides.Add(s);
            }

            stremaReader.Close();
            return peptides;
        }

        public static IEnumerable<string> GetPeptidesFromProtein(string protein, bool fullyTryptic,
                                                              int missedCleavageNumber)
        {
            var peptides = new HashSet<string>();

            for (var sp = 0; sp < protein.Length - 5; sp++)
            {
                var numKR = 0;
                var lterm = sp == 0 || protein[sp - 1] == 'K' || protein[sp - 1] == 'R';

                for (var ep = sp; ep < Math.Min(sp + 50, protein.Length); ep++)
                {
                    var cterm = ep == protein.Length - 1 || protein[ep] == 'K' || protein[ep] == 'R';
                    if (cterm) numKR++;

                    if (ep - sp + 1 >= 5 && numKR - 1 <= missedCleavageNumber &&
                        ((fullyTryptic && lterm && cterm) || !fullyTryptic && (lterm || cterm)))
                    {
                        peptides.Add(protein.Substring(sp, ep - sp + 1));
                    }
                }
            }

            return peptides;
        }

        public static void Shuffle<T>(IList<T> list)
        {
            var rng = new Random();
            var n = list.Count;
            while (n > 1)
            {
                n--;
                var k = rng.Next(n + 1);
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }
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

                            /*var list = new List<char>();
                            foreach (var aa in protein)
                            {
                                list.Add(aa);
                            }
                            Shuffle(list);
                            foreach (var aa in list)
                                revProtein += aa;
                            */
                            for (var i = 0; i < protein.Length; i++)
                            {
                                revProtein += protein[protein.Length - i - 1];
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
