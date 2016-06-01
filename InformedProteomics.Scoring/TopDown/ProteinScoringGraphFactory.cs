using System.Collections.Generic;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.GeneratingFunction;

namespace InformedProteomics.Scoring.TopDown
{
    /*
    public class ProteinScoringGraphFactory
    {

        public ProteinScoringGraphFactory(IMassBinning comparer, AminoAcidSet aminoAcidSet)
        {
            _comparer = comparer;
            _adjList = new List<ScoringGraphEdge>[_comparer.NumberOfBins];

            for (var i = 0; i < _comparer.NumberOfBins; i++) _adjList[i] = new List<ScoringGraphEdge>();

            var terminalModifications = GetTerminalModifications(aminoAcidSet);
            var aminoAcid = GetExtendedAminoAcidSet(aminoAcidSet);

            for (var i = 0; i < _comparer.NumberOfBins; i++)
            {
                var mi = _comparer.GetMass(i);
                var fineNodeMass = mi;

                for (var a = 0; a < aminoAcid.Length; a++)
                {
                    var j = _comparer.GetBinNumber(fineNodeMass + aminoAcid[a].Mass);
                    if (j < 0 || j >= _comparer.NumberOfBins) continue;
                    _adjList[j].Add(new ScoringGraphEdge(i, EdgeScore));
                }

                if (i != 0) continue;
                foreach (var stdAa in AminoAcid.StandardAminoAcidArr)
                {
                    foreach (var terminalMod in terminalModifications)
                    {
                        var modifiedAa = new ModifiedAminoAcid(stdAa, terminalMod);
                        var j = _comparer.GetBinNumber(fineNodeMass + modifiedAa.Mass);
                        if (j < 0 || j >= _comparer.NumberOfBins) continue;
                        _adjList[j].Add(new ScoringGraphEdge(i, EdgeScore));
                    }
                }
            }
        }

        public IScoringGraph CreateScoringGraph(ProductSpectrum deconvSpectrum, double proteinMass)
        {
            if (proteinMass > _comparer.MaxMass || proteinMass < _comparer.MinMass) return null;

            var nodeScores = (deconvSpectrum != null)
                ? GetNodeScores(deconvSpectrum, proteinMass)
                : new int[_comparer.GetBinNumber(proteinMass) + 1];

            var graph = new ProteinScoringGraph(nodeScores, _adjList);
            return graph;
        }

        private static Modification[] GetTerminalModifications(AminoAcidSet aminoAcidSet)
        {
            // N-term, C-term modificaitons
            var terminalModifications = new List<Modification>();
            var residue = AminoAcid.ProteinNTerm.Residue;
            var location = SequenceLocation.ProteinNTerm;
            var aa = aminoAcidSet.GetAminoAcid(residue, location);
            foreach (var modIndex in aminoAcidSet.GetModificationIndices(residue, location))
            {
                var modification = aminoAcidSet.GetModificationParams().GetModification(modIndex);
                terminalModifications.Add(modification);
                //Console.WriteLine(modification.Mass);
            }
            residue = AminoAcid.ProteinCTerm.Residue;
            location = SequenceLocation.ProteinCTerm;
            aa = aminoAcidSet.GetAminoAcid(residue, location);
            foreach (var modIndex in aminoAcidSet.GetModificationIndices(residue, location))
            {
                var modification = aminoAcidSet.GetModificationParams().GetModification(modIndex);
                terminalModifications.Add(modification);
            }

            return terminalModifications.ToArray();
        }

        private static AminoAcid[] GetExtendedAminoAcidSet(AminoAcidSet aaSet)
        {
            var ret = new List<AminoAcid>();
            var modParam = aaSet.GetModificationParams();
            var aminoAcidArray = AminoAcid.StandardAminoAcidArr;

            foreach (var aa in aminoAcidArray)
            {
                ret.Add(aa);
                foreach (var modIndex in aaSet.GetModificationIndices(aa.Residue, SequenceLocation.Everywhere))
                {
                    var aa2 = new ModifiedAminoAcid(aa, modParam.GetModification(modIndex));
                    ret.Add(aa2);
                }
            }

            return ret.ToArray();
        }

        // transform deconvSpectrum to prefix residue map (spectral integer vector)
        private int[] GetNodeScoresByMatchedPeak(ProductSpectrum deconvSpectrum, double proteinMass)
        {
            var baseIonTypes = deconvSpectrum.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            var prefixOffsetMass = baseIonTypes[0].OffsetComposition.Mass;
            var suffixOffsetMass = baseIonTypes[1].OffsetComposition.Mass;

            var numNodes = _comparer.GetBinNumber(proteinMass) + 1;
            var nodeScores = new int[numNodes];

            // assume that peaks are prefixFragment ions
            foreach (var peak in deconvSpectrum.Peaks)
            {
                var prefixIonMass = peak.Mz;
                var prefixFragmentMass = prefixIonMass - prefixOffsetMass;

                var binIndex = _comparer.GetBinNumber(prefixFragmentMass);
                if (binIndex < 0 || binIndex >= numNodes) continue;
                nodeScores[binIndex] = 1;
            }

            // assume that peaks are suffixFragment ions
            foreach (var peak in deconvSpectrum.Peaks)
            {
                var suffixIonMass = peak.Mz;
                var suffixFragmentMass = suffixIonMass - suffixOffsetMass;
                var prefixFragmentMass = proteinMass - suffixFragmentMass;

                var binIndex = _comparer.GetBinNumber(prefixFragmentMass);
                if (binIndex < 0 || binIndex >= numNodes) continue;
                if (nodeScores[binIndex] < 2) nodeScores[binIndex]++;
            }
            return nodeScores;
        }

        private const int EdgeScore = 0;
        private readonly List<ScoringGraphEdge>[] _adjList;
        private readonly IMassBinning _comparer;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static ProteinScoringGraphFactory()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

        internal class AminoAcidInfo
        {
            internal AminoAcidInfo(AminoAcid aa, double prob)
            {
                AminoAcid = aa;
                Probability = prob;
            }
            public AminoAcid AminoAcid { get; private set; }
            public double Probability { get; private set; }
        }

        internal class ProteinScoringGraph : IScoringGraph
        {
            internal ProteinScoringGraph(int[] nodeScores, List<ScoringGraphEdge>[] adjList)
            {
                _nodeScores = nodeScores;
                _adjList = adjList;
            }

            public int GetNodeScore(int nodeIndex)
            {
                return _nodeScores[nodeIndex];
            }

            public IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex)
            {
                return _adjList[nodeIndex];
            }

            public int GetNumNodes()
            {
                return _nodeScores.Length;
            }

            private readonly int[] _nodeScores;
            private readonly List<ScoringGraphEdge>[] _adjList;
        }
    }*/

}
