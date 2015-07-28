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

            _ms2ScanNums = run.GetScanNumbers(2).ToArray();
            _prsmArray = new ProteinSpectrumMatch[_ms2ScanNums.Last() + 1];
            Run = run;
            _massComparerWithBinning = new MzComparerWithBinning(29); // 4ppm
        }

        public void AddIdentificationResult(string path, ProteinSpectrumMatch.SearchTool tool)
        {
            List<ProteinSpectrumMatch> prsmList = null;

            if (tool == ProteinSpectrumMatch.SearchTool.MsAlign)
                prsmList = ReadMsAlignResult(path);
            else if (tool == ProteinSpectrumMatch.SearchTool.MsPathFinder)
                prsmList = ReadMsPathFinderResult(path);

            if (prsmList == null) return;

            foreach (var prsm in prsmList)
            {
                if (_prsmArray[prsm.ScanNum] == null) _prsmArray[prsm.ScanNum] = prsm;
                else
                {
                    if (_prsmArray[prsm.ScanNum].Score < prsm.Score) _prsmArray[prsm.ScanNum] = prsm;
                }
            }
        }
        /*
        public List<ProteinSpectrumMatch> FindByFeature(LcMsFeature feature, Tolerance tolerance)
        {
            return FindByFeature(feature.Mass, feature.MinScanNum, feature.MaxScanNum, tolerance);
        }

        public List<ProteinSpectrumMatch> FindByFeature(double mass, int minScanNum, int maxScanNum, Tolerance tolerance)
        {
            var massTol = tolerance.GetToleranceAsTh(mass);
            var ret = new List<ProteinSpectrumMatch>();

            for (var scanNum = minScanNum + 1; scanNum < maxScanNum; scanNum++)
            {
                List<ProteinSpectrumMatch> prsmList;
                if (_scanToPrSMs.TryGetValue(scanNum, out prsmList))
                {
                    foreach (var prsm in prsmList)
                    {
                        var massDiff = Math.Abs(prsm.Mass - mass);
                        if (massDiff < massTol) ret.Add(prsm);
                    }
                }
            }
            return ret;
        }
        */
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
            var massTol = tolerance.GetToleranceAsTh(mass);
            var massBin = _massComparerWithBinning.GetBinNumber(mass);
            ProteinSpectrumMatch ret = null;
            for (var scanNum = minScanNum + 1; scanNum < maxScanNum; scanNum++)
            {
                if (_prsmArray[scanNum] == null) continue;
                if (Math.Abs(massBin - _massComparerWithBinning.GetBinNumber(_prsmArray[scanNum].Mass)) <= 1)
                {
                    if (ret == null) ret = _prsmArray[scanNum];
                    else
                    {
                        if (ret.Score < _prsmArray[scanNum].Score) ret = _prsmArray[scanNum];
                    }
                }
            }
            return ret;
        }

        public List<ProteinSpectrumMatch> ReadMsAlignResult(string msAlignResultTablePath)
        {
            var parser = new TsvFileParser(msAlignResultTablePath);
            var prsmList = new List<ProteinSpectrumMatch>();

            for (var i = 0; i < parser.NumData; i++)
            {
                var sequence = parser.GetData("Peptide")[i];
                var scanNum = int.Parse(parser.GetData("Scan(s)")[i]);
                var mass = double.Parse(parser.GetData("Adjusted_precursor_mass")[i]);
                var protName = parser.GetData("Protein_name")[i];
                var firstResId = int.Parse(parser.GetData("First_residue")[i]);
                var lastResId = int.Parse(parser.GetData("Last_residue")[i]);
                var score = double.Parse(parser.GetData("#matched_fragment_ions")[i]);
                var sequenceText = parser.GetData("Peptide")[i];

                var fdr = Double.Parse(parser.GetData("FDR")[i]);
                if (fdr > FdrCutoff) continue;

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, protName, firstResId, lastResId, score)
                {
                    SequenceText = sequenceText,
                    SearchToolType = ProteinSpectrumMatch.SearchTool.MsAlign,
                };

                prsmList.Add(prsm);
            }

            return prsmList;
        }

        public List<ProteinSpectrumMatch> ReadMsPathFinderResult(string msPathFinderResultPath)
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

                var prsm = new ProteinSpectrumMatch(sequence, scanNum, mass, protName, firstResId, lastResId, score)
                {
                    SequenceText = sequenceText,
                    SearchToolType = ProteinSpectrumMatch.SearchTool.MsPathFinder,
                };

                prsmList.Add(prsm);

            }

            return prsmList;
        }

        public int CountIdentifiedScans()
        {
            return _prsmArray.Count(prsm => prsm != null);
        }

        public int CountIdentifiedUniqueProteoforms()
        {
            var uniqueSeq = new HashSet<string>();
            foreach (var prsm in _prsmArray.Where(p => p != null))
            {
                uniqueSeq.Add(prsm.SequenceText);
            }
            return uniqueSeq.Count;
        }

        public int CountIdentifiedProteins()
        {
            var uniqueStr = new HashSet<string>();
            foreach (var prsm in _prsmArray.Where(p => p != null))
            {
                uniqueStr.Add(prsm.ProteinName);
            }
            return uniqueStr.Count;
        }

        public readonly LcMsRun Run;
        public int DataId { get; private set; }
        public string DataDesc { get; private set; }
        public const double FdrCutoff = 0.01;
        
        //private readonly Dictionary<int, Dictionary<int, ProteinSpectrumMatch>> _scanMassToPrSM;

        private readonly ProteinSpectrumMatch[] _prsmArray;
        private readonly int[] _ms2ScanNums;

        private readonly MzComparerWithBinning _massComparerWithBinning;


        private static string GetSequenceText(string cleanSequence, string modifications)
        {
            var sequenceText = cleanSequence;
            var parsedModifications = ParseModifications(modifications);
            // Add modifications to sequence
            //parsedModifications.Sort(new CompareModByHighestPosition());   // sort in reverse order for insertion
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
