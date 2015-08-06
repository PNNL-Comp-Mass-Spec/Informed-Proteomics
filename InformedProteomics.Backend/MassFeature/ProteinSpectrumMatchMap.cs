using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Integration;

namespace InformedProteomics.Backend.MassFeature
{
    public class ProteinSpectrumMathMap
    {
        public ProteinSpectrumMathMap(LcMsRun run, int dataid, string dataDesc)
        {
            DataId = dataid;
            DataDesc = dataDesc;
            Run = run;
            _scanNumToMatchMap = new Dictionary<int, ProteinSpectrumMatch>();
        }

        public void LoadIdentificationResult(string path, ProteinSpectrumMatch.SearchTool tool, int maxPrsm = int.MaxValue)
        {
            List<ProteinSpectrumMatch> prsmList = null;

            if (tool == ProteinSpectrumMatch.SearchTool.MsAlign)
                prsmList = ReadMsAlignResult(path, maxPrsm);
            else if (tool == ProteinSpectrumMatch.SearchTool.MsPathFinder)
                prsmList = ReadMsPathFinderResult(path, maxPrsm);

            if (prsmList == null) return;

            ProteinSpectrumMatches = prsmList;
            _scanNumToMatchMap.Clear();

            foreach (var prsm in prsmList)
            {
                if (_scanNumToMatchMap.ContainsKey(prsm.ScanNum))
                {
                    if (_scanNumToMatchMap[prsm.ScanNum].Score < prsm.Score) _scanNumToMatchMap[prsm.ScanNum] = prsm;
                }
                else
                {
                    _scanNumToMatchMap.Add(prsm.ScanNum, prsm);    
                }
            }
        }
        
        public ProteinSpectrumMatch FindByFeature(LcMsFeature feature, Tolerance tolerance)
        {
            return FindByFeature(feature.Mass, feature.MinScanNum, feature.MaxScanNum, tolerance);
        }

        public ProteinSpectrumMatch FindByFeature(double mass, double minTime, double maxTime, Tolerance tolerance)
        {
            var ms1ScanNums = Run.GetMs1ScanVector();
            var minScanNum = -1;
            var maxScanNum = -1;

            for (var i = 1; i < ms1ScanNums.Length; i++)
            {
                var time = Run.GetElutionTime(ms1ScanNums[i]);
                if (minScanNum < 0 && time > minTime)
                {
                    minScanNum = ms1ScanNums[i - 1];
                }
                if (maxScanNum < 0 && time > maxTime)
                {
                    maxScanNum = ms1ScanNums[i];
                    break;
                }
            }
            return FindByFeature(mass, minScanNum, maxScanNum, tolerance);
        }

        public ProteinSpectrumMatch FindByFeature(double mass, int minScanNum, int maxScanNum, Tolerance tolerance)
        {
            var massComparerWithBinning = new MzComparerWithBinning(29); // 4ppm
            var massBin = massComparerWithBinning.GetBinNumber(mass);
            ProteinSpectrumMatch ret = null;
            for (var scanNum = minScanNum + 1; scanNum < maxScanNum; scanNum++)
            {
                if (_scanNumToMatchMap[scanNum] == null) continue;
                if (Math.Abs(massBin - massComparerWithBinning.GetBinNumber(_scanNumToMatchMap[scanNum].Mass)) <= 1)
                {
                    if (ret == null) ret = _scanNumToMatchMap[scanNum];
                    else
                    {
                        if (ret.Score < _scanNumToMatchMap[scanNum].Score) ret = _scanNumToMatchMap[scanNum];
                    }
                }
            }
            return ret;
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
                var protName = parser.GetData("Protein_name")[i];
                var firstResId = int.Parse(parser.GetData("First_residue")[i]);
                var lastResId = int.Parse(parser.GetData("Last_residue")[i]);
                var score = double.Parse(parser.GetData("#matched_fragment_ions")[i]);
                var sequenceText = parser.GetData("Peptide")[i];
                var charge = int.Parse(parser.GetData("Charge")[i]);

                var fdr = Double.Parse(parser.GetData("FDR")[i]);
                if (fdr > FdrCutoff) continue;

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, charge, protName, firstResId, lastResId, score)
                {
                    SequenceText = sequenceText,
                    SearchToolType = ProteinSpectrumMatch.SearchTool.MsAlign,
                };

                prsmList.Add(prsm);

                if (prsmList.Count >= maxPrsm) break;
            }

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsPathFinderResult(string msPathFinderResultPath, int maxPrsm)
        {
            var parser = new TsvFileParser(msPathFinderResultPath);
            var prsmList = new List<ProteinSpectrumMatch>();

            var scoreColumn = parser.GetData("#MatchedFragments") ?? parser.GetData("Score");

            var qValColumn = parser.GetData("QValue");

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Sequence")[i];
                var scanNum = int.Parse(parser.GetData("Scan")[i]);
                var mass = double.Parse(parser.GetData("Mass")[i]);
                var protName = parser.GetData("ProteinName")[i];
                var protDesc = parser.GetData("ProteinDesc")[i];
                var charge = int.Parse(parser.GetData("Charge")[i]);
                //var protDesc = _fastaDb.GetProteinDescription(protName);
                if (protDesc != null) protName = protName + " " + protDesc;

                var firstResId = int.Parse(parser.GetData("Start")[i]);
                var lastResId = int.Parse(parser.GetData("End")[i]);
                var score = int.Parse(scoreColumn[i]);
                var mod = parser.GetData("Modifications")[i];

                if (qValColumn != null)
                {
                    var fdr = double.Parse(qValColumn[i]);
                    if (fdr > FdrCutoff) continue;
                }
                
                //var seqKeyPair = SetModifications(sequence, mod);
                var sequenceText = GetSequenceText(sequence, mod);

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, charge, protName, firstResId, lastResId, score)
                {
                    SequenceText = sequenceText,
                    SearchToolType = ProteinSpectrumMatch.SearchTool.MsPathFinder,
                };

                prsmList.Add(prsm);
                
                if (prsmList.Count >= maxPrsm) break;
            }

            return prsmList;
        }

        public int CountIdentifiedScans()
        {
            return _scanNumToMatchMap.Count;
        }

        public int CountIdentifiedUniqueProteoforms()
        {
            if (ProteinSpectrumMatches == null) return 0;
            var uniqueSeq = new HashSet<string>();
            foreach (var prsm in ProteinSpectrumMatches) uniqueSeq.Add(prsm.SequenceText);
            return uniqueSeq.Count;
        }

        public int CountIdentifiedProteins()
        {
            if (ProteinSpectrumMatches == null) return 0;
            var uniqueSeq = new HashSet<string>();
            foreach (var prsm in ProteinSpectrumMatches) uniqueSeq.Add(prsm.ProteinName);
            return uniqueSeq.Count;
        }

        public List<ProteinSpectrumMatch> ProteinSpectrumMatches { get; private set; }
        
        public readonly LcMsRun Run;
        public int DataId { get; private set; }
        public string DataDesc { get; private set; }
        public const double FdrCutoff = 0.01;
        private readonly Dictionary<int, ProteinSpectrumMatch> _scanNumToMatchMap;
        
        private static string GetSequenceText(string cleanSequence, string modifications)
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

        private static List<Tuple<int, string>> ParseModifications(string modifications)
        {
            var mods = modifications.Split(',');
            var parsedMods = new List<Tuple<int, string>>();
            if (mods.Length < 1 || mods[0] == string.Empty)
            {
                return parsedMods;
            }

            foreach (var modParts in mods.Select(mod => mod.Split(' ')))
            {
                if (modParts.Length < 0)
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
