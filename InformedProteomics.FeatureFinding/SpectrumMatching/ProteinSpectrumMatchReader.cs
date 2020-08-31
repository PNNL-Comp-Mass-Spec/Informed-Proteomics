using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.SearchResults;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.FeatureFinding.SpectrumMatching
{
    public class ProteinSpectrumMatchReader
    {
        // ReSharper disable once CommentTypo
        // Ignore Spelling: Desc, msgfdb

        public readonly double FdrCutoff;

        public ProteinSpectrumMatchReader(double fdrCutoff = 0.01)
        {
            FdrCutoff = fdrCutoff;
        }

        public List<ProteinSpectrumMatch> LoadIdentificationResult(string path, ProteinSpectrumMatch.SearchTool tool = ProteinSpectrumMatch.SearchTool.Unknown, int maxPrsm = int.MaxValue)
        {
            List<ProteinSpectrumMatch> prsmList = null;

            if (tool == ProteinSpectrumMatch.SearchTool.Unknown)
            {
                if (path.EndsWith("IcTda.tsv", StringComparison.OrdinalIgnoreCase))
                    tool = ProteinSpectrumMatch.SearchTool.MsPathFinder;
                else if (path.EndsWith("MSAlign_ResultTable.txt", StringComparison.OrdinalIgnoreCase))
                    tool = ProteinSpectrumMatch.SearchTool.MsAlign;
                // ReSharper disable StringLiteralTypo
                else if (path.EndsWith("msgfplus_syn.txt", StringComparison.OrdinalIgnoreCase) ||
                         path.EndsWith("msgfplus_fht.txt", StringComparison.OrdinalIgnoreCase) ||
                         path.EndsWith("msgfdb_syn.txt", StringComparison.OrdinalIgnoreCase) ||
                         path.EndsWith("msgfdb_fht.txt", StringComparison.OrdinalIgnoreCase))
                    tool = ProteinSpectrumMatch.SearchTool.MsGfPlus;
                // ReSharper restore StringLiteralTypo
            }

            if (tool == ProteinSpectrumMatch.SearchTool.MsAlign)
                prsmList = ReadMsAlignResult(path, maxPrsm);
            else if (tool == ProteinSpectrumMatch.SearchTool.MsPathFinder)
                prsmList = ReadMsPathFinderResult(path, maxPrsm);
            else if (tool == ProteinSpectrumMatch.SearchTool.MsGfPlus)
                prsmList = ReadMsGfPlusResult(path, maxPrsm);

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsAlignResult(string msAlignResultTablePath, int maxPrsm)
        {
            var parser = new TsvFileParser(msAlignResultTablePath);
            var prsmList = new List<ProteinSpectrumMatch>();

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Peptide")[i];
                var scanNum = int.Parse(parser.GetData("Scan(s)")[i]);
                var mass = double.Parse(parser.GetData("Precursor_mass")[i]);
                var protNameDesc = parser.GetData("Protein_name")[i];

                var k = protNameDesc.IndexOf(' ');
                var protName = (k < 0) ? protNameDesc : protNameDesc.Substring(0, k);
                var protDesc = (k < 0) ? protNameDesc : protNameDesc.Substring(k+1);

                var firstResId = int.Parse(parser.GetData("First_residue")[i]);
                var lastResId = int.Parse(parser.GetData("Last_residue")[i]);
                var score = double.Parse(parser.GetData("#matched_fragment_ions")[i]);
                var sequenceText = parser.GetData("Peptide")[i];
                var charge = int.Parse(parser.GetData("Charge")[i]);
                var eValue = double.Parse(parser.GetData("E-value")[i]);

                var fdr = double.Parse(parser.GetData("FDR")[i]);
                if (fdr > FdrCutoff) continue;

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, charge, protName, protDesc, firstResId, lastResId, score, ProteinSpectrumMatch.SearchTool.MsAlign)
                {
                    SequenceText = sequenceText,
                    SpectralEValue = eValue,
                };

                prsmList.Add(prsm);

                if (prsmList.Count >= maxPrsm) break;
            }

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsGfPlusResult(string msgfResultPath, int maxPrsm)
        {
            var parser = new TsvFileParser(msgfResultPath);
            var prsmList = new List<ProteinSpectrumMatch>();
            var prevScanNum = -1;

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Peptide")[i];
                var scanNum = int.Parse(parser.GetData("Scan")[i]);

                if (prevScanNum == scanNum) continue;
                prevScanNum = scanNum;

                var mz = double.Parse(parser.GetData("PrecursorMZ")[i]);
                var protName = parser.GetData("Protein")[i];
                var protDesc = "";
                var score = double.Parse(parser.GetData("MSGFScore")[i]);
                var charge = int.Parse(parser.GetData("Charge")[i]);

                var seq = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequence);
                var sequenceText = GetSequenceText(seq);
                var mass = (mz - Constants.Proton)*charge;
                var firstResId = 0;
                var lastResId = 0;
                var fdr = double.Parse(parser.GetData("QValue")[i]);
                if (fdr > FdrCutoff) continue;

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, charge, protName, protDesc, firstResId, lastResId, score, ProteinSpectrumMatch.SearchTool.MsGfPlus)
                {
                    SequenceText = sequenceText,
                };

                prsmList.Add(prsm);

                if (prsmList.Count >= maxPrsm) break;
            }

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsPathFinderResultOld(string msPathFinderResultPath, int maxPrsm, double minScore = 3, double maxScore = int.MaxValue)
        {
            var parser = new TsvFileParser(msPathFinderResultPath);
            var prsmList = new List<ProteinSpectrumMatch>();

            var scoreColumn = parser.GetData("#MatchedFragments") ?? parser.GetData("Score");
            var qValColumn = parser.GetData("QValue");

            var eValueColumn = parser.GetData("SpecEValue");

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Sequence")[i];
                var scanNum = int.Parse(parser.GetData("Scan")[i]);
                var mass = double.Parse(parser.GetData("Mass")[i]);
                var protName = parser.GetData("ProteinName")[i];
                var protDesc = parser.GetData("ProteinDesc")[i];
                var charge = int.Parse(parser.GetData("Charge")[i]);

                var firstResId = int.Parse(parser.GetData("Start")[i]);
                var lastResId = int.Parse(parser.GetData("End")[i]);
                var score = double.Parse(scoreColumn[i]);
                var mod = parser.GetData("Modifications")[i];
                var eValue = (eValueColumn != null) ? double.Parse(parser.GetData("SpecEValue")[i]) : 0;

                var pre = parser.GetData("Pre")[i];
                var post = parser.GetData("Post")[i];
                var proteinLen = int.Parse(parser.GetData("ProteinLength")[i]);

                if (score < minScore || score > maxScore) continue;

                if (qValColumn != null)
                {
                    var fdr = double.Parse(qValColumn[i]);
                    if (fdr > FdrCutoff) continue;
                }

                var sequenceText = GetSequenceText(sequence, mod);

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, charge, protName, protDesc, firstResId, lastResId, score, ProteinSpectrumMatch.SearchTool.MsPathFinder)
                {
                    SequenceText = sequenceText,
                    Modifications = mod,
                    Pre = pre,
                    Post = post,
                    ProteinLength = proteinLen,
                    SpectralEValue = eValue,
                };

                prsmList.Add(prsm);

                if (prsmList.Count >= maxPrsm) break;
            }

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsPathFinderResult(string msPathFinderResultPath, int maxPrsm, double minScore = 3, double maxScore = int.MaxValue)
        {
            var prsmList = new List<ProteinSpectrumMatch>();

            foreach (var result in DatabaseSearchResultData.ReadResultsFromFile(msPathFinderResultPath))
            {
                if (result.NumMatchedFragments < minScore || maxScore < result.NumMatchedFragments)
                {
                    continue;
                }

                if (result.HasTdaScores && result.QValue > FdrCutoff)
                {
                    continue;
                }

                var prsm = new ProteinSpectrumMatch(result.Sequence, result.ScanNum, result.Mass, result.Charge, result.ProteinName, result.ProteinDescription, result.Start, result.End, result.NumMatchedFragments, ProteinSpectrumMatch.SearchTool.MsPathFinder)
                {
                    SequenceText = GetSequenceText(result.Sequence, result.Modifications),
                    Modifications = result.Modifications,
                    Pre = result.Pre,
                    Post = result.Post,
                    ProteinLength = result.ProteinLength,
                    SpectralEValue = result.SpecEValue,
                };

                prsmList.Add(prsm);

                if (prsmList.Count >= maxPrsm)
                    break;
            }

            return prsmList;
        }

        public static string GetSequenceText(Sequence sequence)
        {
            var sequenceText = new StringBuilder();
            foreach (var aa in sequence)
            {
                sequenceText.Append(aa.Residue);

                if (aa is ModifiedAminoAcid ma)
                {
                    var modLabel = string.Format("[{0}]", ma.Modification.Name);
                    sequenceText.Append(modLabel);
                }
            }
            return sequenceText.ToString();
        }

        public static string GetSequenceText(string cleanSequence, string modifications)
        {
            var sequenceText = cleanSequence;
            var parsedModifications = ParseModifications(modifications);
            foreach (var mod in parsedModifications.OrderByDescending(m => m.Item1))
            {
                var modLabel = string.Format("[{0}]", mod.Item2);
                sequenceText = sequenceText.Insert(mod.Item1, modLabel);
            }
            return sequenceText;
        }

        private static IEnumerable<Tuple<int, string>> ParseModifications(string modifications)
        {
            var mods = modifications.Split(',');
            var parsedMods = new List<Tuple<int, string>>();
            if (mods.Length < 1 || mods[0] == string.Empty)
            {
                return parsedMods;
            }

            foreach (var modParts in mods.Select(mod => mod.Split(' ')))
            {
                if (modParts.Length == 0)
                {
                    throw new FormatException("Unknown Modification");
                }

                var modName = modParts[0];
                var modPos = Convert.ToInt32(modParts[1]);
                parsedMods.Add(new Tuple<int, string>(modPos, modName));
            }
            return parsedMods;
        }
    }
}
