using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class PrSmContainer
    {
        public PrSmContainer(int dataid, string dataDesc, string fastaFilePath)
        {
            DataId = dataid;
            DataDesc = dataDesc;

            _fastaDb = new FastaDatabase(fastaFilePath);
            _fastaDb.Read();
            _scanToPrSMs = new Dictionary<int, List<ProteinSpectrumMatch>>();
            Fdr1Threshold = double.MaxValue;
        }

        public void AddIdentificationResult(string path, ProteinSpectrumMatch.SearchTool tool)
        {
            List<ProteinSpectrumMatch> prsmList = null;

            if (tool == ProteinSpectrumMatch.SearchTool.MsAlign)
                prsmList = ReadMsAlignResult(path);
            else if (tool == ProteinSpectrumMatch.SearchTool.MsPathFinder)
                prsmList = ReadMsPathFinderResult(path);

            if (prsmList != null)
            {
                MergeIdentifcationResult(prsmList);

                var uniqueSeq = new HashSet<string>();
                foreach (var prsm in prsmList)
                {
                    uniqueSeq.Add(prsm.SequenceText);
                }
                //Console.WriteLine("PrSMs = {0}; unique proteins = {1}", prsmList.Count, uniqueSeq.Count);
            }
        }

        public List<ProteinSpectrumMatch> FindByFeature(LcMsFeature feature, Tolerance tolerance)
        {
            var massTol = tolerance.GetToleranceAsTh(feature.Mass);
            var ret = new List<ProteinSpectrumMatch>();

            for (var scanNum = feature.MinScanNum + 1; scanNum < feature.MaxScanNum; scanNum++)
            {
                List<ProteinSpectrumMatch> prsmList;
                if (_scanToPrSMs.TryGetValue(scanNum, out prsmList))
                {
                    foreach (var prsm in prsmList)
                    {
                        var massDiff = Math.Abs(prsm.Mass - feature.Mass);
                        if (massDiff < massTol) ret.Add(prsm);
                    }
                }
            }

            return ret;
        }

        private void MergeIdentifcationResult(IEnumerable<ProteinSpectrumMatch> newPrsmList)
        {
            foreach (var prsm in newPrsmList)
            {
                List<ProteinSpectrumMatch> prsmList;
                if (_scanToPrSMs.TryGetValue(prsm.ScanNum, out prsmList))
                {
                    if (!prsmList.Contains(prsm)) prsmList.Add(prsm);
                }
                else
                {
                    _scanToPrSMs.Add(prsm.ScanNum, new List<ProteinSpectrumMatch> {prsm});
                }
            }
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
                if (fdr > 0.01) continue;

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
                var protDesc = _fastaDb.GetProteinDescription(protName);
                if (protDesc != null) protName = protName + " " + protDesc;

                var firstResId = int.Parse(parser.GetData("Start")[i]);
                var lastResId = int.Parse(parser.GetData("End")[i]);
                var score = double.Parse(scoreColumn[i]);
                var mod = parser.GetData("Modifications")[i];

                if (qValColumn != null)
                {
                    var fdr = double.Parse(qValColumn[i]);
                    if (fdr > 0.01) continue;
                    if (score < Fdr1Threshold) Fdr1Threshold = score;
                }
                else
                {
                    if (score < Fdr1Threshold) continue;
                }

                var seqKeyPair = SetModifications(sequence, mod);
                var sequenceText = seqKeyPair.Item2;

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
            return _scanToPrSMs.Count;
        }

        public int CountIdentifiedUniqueProteoforms()
        {
            var uniqueSeq = new HashSet<string>();
            foreach (var prsmList in _scanToPrSMs.Values)
            {
                foreach(var prsm in prsmList) uniqueSeq.Add(prsm.SequenceText);
            }
            return uniqueSeq.Count;
        }

        public int CountIdentifiedProteins()
        {
            return _scanToPrSMs.Values.Sum(prsmList => prsmList.Count);
        }

        private static Tuple<Sequence, string> SetModifications(string cleanSequence, string modifications)
        {
            // Build Sequence AminoAcid list
            var sequence = new Sequence(cleanSequence, new AminoAcidSet());
            var sequenceText = cleanSequence;
            var parsedModifications = ParseModifications(modifications);

            // Add modifications to sequence
            parsedModifications.Sort(new CompareModByHighestPosition());   // sort in reverse order for insertion
            foreach (var mod in parsedModifications)
            {
                var pos = mod.Item1;
                if (pos > 0)
                {
                    pos--;
                }

                var modLabel = string.Format("[{0}]", mod.Item2.Name);
                sequenceText = sequenceText.Insert(mod.Item1, modLabel);
                var aa = sequence[pos];
                var modaa = new ModifiedAminoAcid(aa, mod.Item2);
                sequence[pos] = modaa;
            }

            return new Tuple<Sequence, string>(new Sequence(sequence), sequenceText);
        }

        private static List<Tuple<int, Modification>> ParseModifications(string modifications)
        {
            var mods = modifications.Split(',');
            var parsedMods = new List<Tuple<int, Modification>>();
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
                var modification = Modification.Get(modName);
                parsedMods.Add(new Tuple<int, Modification>(modPos, modification));
                if (modification == null)
                {
                    throw new Exception(string.Format("Found an unrecognized modification: {0}", modName));
                }
            }

            return parsedMods;
        }

        public int DataId { get; private set; }
        public string DataDesc { get; private set; }

        public double Fdr1Threshold { get; private set; }
        private readonly Dictionary<int, List<ProteinSpectrumMatch>> _scanToPrSMs;
        private readonly FastaDatabase _fastaDb;

        internal class CompareModByHighestPosition : IComparer<Tuple<int, Modification>>
        {
            public int Compare(Tuple<int, Modification> x, Tuple<int, Modification> y)
            {
                return y.Item1.CompareTo(x.Item1);
            }
        }    

    }
}
