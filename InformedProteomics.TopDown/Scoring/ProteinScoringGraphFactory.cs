using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    public class ProteinScoringGraphFactory
    {
        public ProteinScoringGraphFactory(IMassBinning comparer, AminoAcidSet aminoAcidSet)
        {
            _comparer = comparer;
            _adjList = new LinkedList<ScoringGraphEdge>[_comparer.NumberOfBins];

            for (var i = 0; i < _comparer.NumberOfBins; i++)
            {
                _adjList[i] = new LinkedList<ScoringGraphEdge>();
            }

            var terminalModifications = FilteredProteinMassBinning.GetTerminalModifications(aminoAcidSet);
            var aminoAcidArray = FilteredProteinMassBinning.GetExtendedAminoAcidArray(aminoAcidSet);

            for (var i = 0; i < _comparer.NumberOfBins; i++)
            {
                var mi = _comparer.GetMass(i);
                var fineNodeMass = mi;

                foreach (var aa in aminoAcidArray)
                {
                    var j = _comparer.GetBinNumber(fineNodeMass + aa.Mass);
                    if (j < 0 || j >= _comparer.NumberOfBins)
                    {
                        continue;
                    }

                    _adjList[j].AddLast(new ScoringGraphEdge(i));

                    if (i == 0 && aa is not ModifiedAminoAcid)
                    {
                        foreach (var terminalMod in terminalModifications)
                        {
                            var modifiedAa = new ModifiedAminoAcid(aa, terminalMod);
                            j = _comparer.GetBinNumber(fineNodeMass + modifiedAa.Mass);
                            if (j < 0 || j >= _comparer.NumberOfBins)
                            {
                                continue;
                            }

                            _adjList[j].AddLast(new ScoringGraphEdge(i));
                        }
                    }
                }
            }
        }

        public IScoringGraph CreateScoringGraph(CompositeScorerBasedOnDeconvolutedSpectrum scorer, double proteinMass)
        {
            if (proteinMass > _comparer.MaxMass || proteinMass < _comparer.MinMass)
            {
                return null;
            }

            var nodeScores = scorer.GetNodeScores(proteinMass);
            var graph = new ProteinScoringGraph(nodeScores[0], nodeScores[1], _adjList, _comparer);

            return graph;
        }
        private readonly LinkedList<ScoringGraphEdge>[] _adjList;
        private readonly IMassBinning _comparer;

        //private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        // static ProteinScoringGraphFactory()
        // {
        //     BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
        //     BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        // }

        internal class ProteinScoringGraph : IScoringGraph
        {
            internal ProteinScoringGraph(double?[] nodeScoresByPrefixIon, double?[] nodeScoresBySuffixIon, LinkedList<ScoringGraphEdge>[] adjList, IMassBinning comparer)
            {
                _nodeScoresByPrefixIon = nodeScoresByPrefixIon;
                _nodeScoresBySuffixIon = nodeScoresBySuffixIon;
                _adjList = adjList;
                _comparer = comparer;
            }

            public double GetNodeScore(int nodeIndex)
            {
                return (_nodeScoresByPrefixIon[nodeIndex] ?? 0d) + (_nodeScoresBySuffixIon[nodeIndex] ?? 0d);
            }

            public double GetEdgeScore(int nodeIndex1, int nodeIndex2)
            {
                var edgeScore = 0d;
                if (_nodeScoresByPrefixIon[nodeIndex1] != null && _nodeScoresByPrefixIon[nodeIndex2] != null)
                {
                    edgeScore += CompositeScorer.ScoreParam.Prefix.ConsecutiveMatch;
                }
                if (_nodeScoresBySuffixIon[nodeIndex1] != null && _nodeScoresBySuffixIon[nodeIndex2] != null)
                {
                    edgeScore += CompositeScorer.ScoreParam.Suffix.ConsecutiveMatch;
                }
                return edgeScore;
            }

            public IEnumerable<IScoringGraphEdge> GetEdges(int nodeIndex)
            {
                return nodeIndex >= GetNumNodes() ? Enumerable.Empty<ScoringGraphEdge>() : _adjList[nodeIndex];
            }

            public int GetNumNodes()
            {
                return _nodeScoresByPrefixIon.Length;
            }

            public double ScoreSequence(Sequence sequence)
            {
                var cleavages = sequence.GetInternalCleavages();
                var score = 0.0;
                var prevNTermBin = 0;
                foreach (var cleavage in cleavages)
                {
                    var nTermBin = _comparer.GetBinNumber(cleavage.PrefixComposition.Mass);

                    if (nTermBin < 0)
                    {
                        prevNTermBin = nTermBin;
                        continue;
                    }

                    score += GetNodeScore(nTermBin);
                    score += prevNTermBin >= 0 ? GetEdgeScore(prevNTermBin, nTermBin) : 0;

                    prevNTermBin = nTermBin;
                }

                return score;
            }

            private readonly IMassBinning _comparer;
            private readonly double?[] _nodeScoresByPrefixIon;
            private readonly double?[] _nodeScoresBySuffixIon;
            private readonly LinkedList<ScoringGraphEdge>[] _adjList;
        }
        /*
        // transform deconvSpectrum to prefix residue map (spectral integer vector)
        private int[] GetNodeScoresByMatchedPeak(DeconvolutedSpectrum deconvSpectrum, double proteinMass)
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
        }*/
    }
}
