using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using PRISM;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Creates a graph for scoring a sequence or annotation
    /// </summary>
    public class SequenceGraph
    {
        /// <summary>
        /// Create a graph representing the annotation. Annotation is reversed.
        /// </summary>
        /// <param name="aaSet">amino acid set</param>
        /// <param name="annotation">annotation (e.g. G.PEPTIDER.K or _.PEPTIDER._)</param>
        /// <returns></returns>
        public static SequenceGraph CreateGraph(AminoAcidSet aaSet, string annotation)
        {
            const char delimiter = (char)FastaDatabaseConstants.Delimiter;
            if (annotation == null || !Regex.IsMatch(annotation, @"^[A-Z" + delimiter + @"]\.[A-Z]+\.[A-Z" + delimiter + @"]$")) return null;

            var nTerm = annotation[0] == FastaDatabaseConstants.Delimiter
                                  ? AminoAcid.ProteinNTerm
                                  : AminoAcid.PeptideNTerm;
            var cTerm = annotation[annotation.Length - 1] == FastaDatabaseConstants.Delimiter
                                  ? AminoAcid.ProteinCTerm
                                  : AminoAcid.PeptideCTerm;

            var sequence = annotation.Substring(2, annotation.Length - 4);
            return CreateGraph(aaSet, nTerm, sequence, cTerm);
        }

        /// <summary>
        /// Create a graph representing the sequence. Sequence is reversed.
        /// </summary>
        /// <param name="aaSet">amino acid set</param>
        /// <param name="nTerm">N-term amino acid</param>
        /// <param name="sequence">sequence</param>
        /// <param name="cTerm">C-term amino acid</param>
        /// <returns>sequence graph</returns>
        public static SequenceGraph CreateGraph(AminoAcidSet aaSet, AminoAcid nTerm, string sequence, AminoAcid cTerm)
        {
            var seqGraph = new SequenceGraph(aaSet, nTerm, sequence, cTerm);

            return seqGraph.IsValid ? seqGraph : null;
        }

        /// <summary>
        /// Amino acid set used in the graph
        /// </summary>
        public AminoAcidSet AminoAcidSet { get; }

        /// <summary>
        /// Modifications used in this graph
        /// </summary>
        public ModificationParams ModificationParams => _modificationParams;

        /// <summary>
        /// If the sequence graph is valid
        /// </summary>
        public bool IsValid { get; }

        /// <summary>
        /// Number of N-Terminus cleavages
        /// </summary>
        public int NumNTermCleavages { get; private set; }

        /// <summary>
        /// Gets the number of possible compositions of the current sequence
        /// </summary>
        /// <returns>the number of possible compositions</returns>
        public int GetNumProteoformCompositions()
        {
            return _graph[_index].Length;
        }

        /// <summary>
        /// Gets the number of distinct compositions of the current sequence
        /// </summary>
        /// <returns>the number of possible compositions</returns>
        public int GetNumDistinctSequenceCompositions()
        {
            var compositions = new HashSet<Composition.Composition>();
            for (var nodeIndex = 0; nodeIndex < _graph[_index].Length; nodeIndex++)
            {
                compositions.Add(GetComposition(_index, nodeIndex));
            }
            return compositions.Count;
        }

        /// <summary>
        /// Gets all possible compositions of the current sequence
        /// </summary>
        /// <returns>all possible compositions</returns>
        public Composition.Composition[] GetSequenceCompositions()
        {
            var numCompositions = _graph[_index].Length;
            var compositions = new Composition.Composition[numCompositions];
            for (var modIndex = 0; modIndex < numCompositions; modIndex++)
            {
                compositions[modIndex] = GetComposition(_index, modIndex);
            }
            return compositions;
        }

        /// <summary>
        /// Gets number of possible proteoforms
        /// </summary>
        /// <returns>number of possible proteoforms</returns>
        public int GetNumProteoforms()
        {
            return _graph[_index].Length;
        }

        /// <summary>
        /// Get the number of proteoform sequences possible with the provided number of dynamic modifications
        /// </summary>
        /// <param name="nDynMods"></param>
        /// <returns></returns>
        public double GetNumProteoformSequencesByNumMods(int nDynMods)
        {
            if (nDynMods < 1 || nDynMods > ModificationParams.MaxNumDynModsPerSequence) return 0;

            var numModCombs = _graph[_index].Length;
            var nProteoforms = 0d;
            for (var modIndex = 0; modIndex < numModCombs; modIndex++)
            {
                var modCombs = ModificationParams.GetModificationCombination(_graph[_index][modIndex].ModificationCombinationIndex);
                if (modCombs.GetNumModifications() < nDynMods) continue;
                if (modCombs.GetNumModifications() > nDynMods) break;

                nProteoforms += GetNumProteoformSequences(modIndex);
                //if (nProteoforms > double.MaxValue) return double.MaxValue; // overflow
            }
            return nProteoforms;
        }

        /// <summary>
        /// Get number of possible proteoform sequences for the speicified modification combination
        /// </summary>
        /// <param name="modIndex">index of modification combination</param>
        /// <returns>number of sequences</returns>
        public double GetNumProteoformSequences(int modIndex)
        {
            var numModCombs = _graph[_index].Length;
            if (modIndex < 0 || modIndex >= numModCombs) return 0;

            if (_countTable == null)
            {
                _countTable = new double[_maxSeqIndex][];
                for (var i = 0; i < _maxSeqIndex; i++)
                {
                    _countTable[i] = new double[_graph[i].Length];
                }

                _countTable[0][0] = 1;
                for (var i = 1; i < _maxSeqIndex; i++)
                {
                    for (var j = 0; j < _graph[i].Length; j++)
                    {
                        foreach (var k in _graph[i][j].GetPrevNodeIndices())
                        {
                            _countTable[i][j] += _countTable[i - 1][k];
                        }
                    }
                }
            }

            return _countTable[_index][modIndex];
        }

        private double[][] _countTable;

        /// <summary>
        /// Perform  N-Terminus cleavage on the sequence graph
        /// </summary>
        public void CleaveNTerm()
        {
            _index = _index - 3;
            for (var seqIndex = _index + 1; seqIndex < _index + 3; seqIndex++)
            {
                for (var modIndex = 0; modIndex < _graph[seqIndex].Length; modIndex++)
                {
                    _nodeComposition[seqIndex][modIndex] = null;
                }
            }

            ++NumNTermCleavages;
            SetNTerminalAminoAcid(_nTerm);
            AddAminoAcid(_sequence[NumNTermCleavages]);
            AddAminoAcid(_nTerm.Residue);
        }

        /// <summary>
        /// Get the modification combinations in the current sequence graph
        /// </summary>
        /// <returns></returns>
        public ModificationCombination[] GetModificationCombinations()
        {
            var numModCombs = _graph[_index].Length;
            var modCombs = new ModificationCombination[numModCombs];
            for (var modIndex = 0; modIndex < numModCombs; modIndex++)
            {
                modCombs[modIndex] = ModificationParams.GetModificationCombination(_graph[_index][modIndex].ModificationCombinationIndex);
            }
            return modCombs;
        }

        /// <summary>
        /// Gets the composition of the sequence without variable modification.
        /// </summary>
        /// <returns></returns>
        public Composition.Composition GetUnmodifiedSequenceComposition()
        {
            return _suffixComposition[_index];
        }

        /// <summary>
        /// Gets the number of possible product compositions of the current sequence
        /// </summary>
        /// <returns>the number of possible product compositions</returns>
        public int GetNumFragmentCompositions()
        {
            var numNodes = 0;
            for (var index = 2; index < _index - 2; index++)
            {
                numNodes += _graph[index].Length;
            }
            return numNodes;
        }

        /// <summary>
        /// Get the compositions of all fragment nodes
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Composition.Composition> GetAllFragmentNodeCompositions()
        {
            for (var seqIndex = 2; seqIndex < _maxSeqIndex - 2; seqIndex++)
            {
                for (var modIndex = 0; modIndex < _graph[seqIndex].Length; modIndex++)
                {
                    yield return GetComposition(seqIndex, modIndex);
                }
            }
        }

        /// <summary>
        /// Get the fragment compositions for the provided index and number of N-Terminus cleavages
        /// </summary>
        /// <param name="modIndex"></param>
        /// <param name="numNTermCleavages"></param>
        /// <returns></returns>
        public IEnumerable<Composition.Composition> GetFragmentCompositions(int modIndex, int numNTermCleavages)
        {
            var seqIndex = _index - 1 - numNTermCleavages;

            var fragmentCompositions = new List<Composition.Composition>();

            var modIndexList = new HashSet<int> { modIndex };
            for (var si = seqIndex; si >= 2; si--)
            {
                var newModIndexList = new HashSet<int>();
                foreach (var mi in modIndexList)
                {
                    if (si < seqIndex) fragmentCompositions.Add(GetComposition(si, mi));
                    var node = _graph[si][mi];
                    foreach (var prevModIndex in node.GetPrevNodeIndices())
                    {
                        newModIndexList.Add(prevModIndex);
                    }
                }
                modIndexList = newModIndexList;
            }

            return fragmentCompositions;
        }

        /// <summary>
        /// Set the provided index as the sink
        /// </summary>
        /// <param name="modIndex"></param>
        public void SetSink(int modIndex)
        {
            _sinkModIndex = modIndex;
            _sinkSequenceComposition = GetComposition(_index, modIndex);
            _sinkSequenceCompositionWithH2O = _sinkSequenceComposition + Composition.Composition.H2O;
            foreach (var cn in _compNodeComposition) Array.Clear(cn, 0, cn.Length);
        }

        /// <summary>
        /// Get the Sink sequence composition with H2O added
        /// </summary>
        /// <returns></returns>
        public Composition.Composition GetSinkSequenceCompositionWithH2O()
        {
            return _sinkSequenceCompositionWithH2O;
        }

        private Composition.Composition _sinkSequenceComposition;
        private Composition.Composition _sinkSequenceCompositionWithH2O;
        private int _sinkModIndex;
        private readonly Composition.Composition[][] _compNodeComposition;

        /// <summary>
        /// Get the fragment score using the provided scorer
        /// </summary>
        /// <param name="scorer"></param>
        /// <returns></returns>
        public double GetFragmentScore(IScorer scorer)
        {
            var nodeScore = new double?[_maxSeqIndex][];
            var maxScore = new double?[_maxSeqIndex][];
            for (var si = 0; si < _maxSeqIndex; si++)
            {
                nodeScore[si] = new double?[_graph[si].Length];
                maxScore[si] = new double?[_graph[si].Length];
            }
            maxScore[0][0] = 0.0;
            nodeScore[_index][_sinkModIndex] = 0.0;

            var fragmentScore = GetFragmentScore(_index, _sinkModIndex, scorer, nodeScore, maxScore);
            return fragmentScore;
        }

        /// <summary>
        /// Get the fragment score and string-format modifications of this sequence graph
        /// </summary>
        /// <param name="scorer"></param>
        /// <returns></returns>
        public Tuple<double, string> GetFragmentScoreAndModifications(IScorer scorer)
        {
            var nodeScore = new double?[_maxSeqIndex][];
            var maxScore = new Tuple<double, string>[_maxSeqIndex][];
            for (var si = 0; si < _maxSeqIndex; si++)
            {
                nodeScore[si] = new double?[_graph[si].Length];
                maxScore[si] = new Tuple<double, string>[_graph[si].Length];
            }
            maxScore[0][0] = new Tuple<double, string>(0.0, "");
            nodeScore[_index][_sinkModIndex] = 0.0;

            var fragmentScore = GetFragmentScoreAndModifications(_index, _sinkModIndex, scorer, nodeScore, maxScore);

            return new Tuple<double, string>(fragmentScore.Item1, fragmentScore.Item2);
        }

        /// <summary>
        /// Generate a sequence graph for the provided data
        /// </summary>
        /// <param name="aminoAcidSet"></param>
        /// <param name="nTerm"></param>
        /// <param name="sequence"></param>
        /// <param name="cTerm"></param>
        protected SequenceGraph(AminoAcidSet aminoAcidSet, AminoAcid nTerm, string sequence, AminoAcid cTerm)
        {
            AminoAcidSet = aminoAcidSet;
            _sequence = sequence;
            _nTerm = nTerm;

            _modificationParams = aminoAcidSet.GetModificationParams();

            _maxSeqIndex = sequence.Length + 3;  // init + C-term + sequence length + N-term

            _index = 0;

            _aminoAcidSequence = new AminoAcid[_maxSeqIndex];
            _aminoAcidSequence[0] = AminoAcid.Empty;

            _suffixComposition = new Composition.Composition[_maxSeqIndex];
            _suffixComposition[0] = Composition.Composition.Zero;

            _graph = new Node[_maxSeqIndex][];
            _graph[0] = new[] { new Node(0) };

            _nodeComposition = new Composition.Composition[_maxSeqIndex][]; //, _modificationParams.NumModificationCombinations];
            _compNodeComposition = new Composition.Composition[_maxSeqIndex][]; //, _modificationParams.NumModificationCombinations];

            for (var i = 0; i < _maxSeqIndex; i++)
            {
                _compNodeComposition[i] = new Composition.Composition[_modificationParams.NumModificationCombinations];
                _nodeComposition[i] = new Composition.Composition[_modificationParams.NumModificationCombinations];
            }

            NumNTermCleavages = 0;
            IsValid = true;
            SetNTerminalAminoAcid(nTerm);
            AddAminoAcid(cTerm.Residue);
            for (var i = sequence.Length - 1; i >= 0; i--)
            {
                if (AddAminoAcid(sequence[i]) == false)
                {
                    IsValid = false;
                    break;
                }
            }
            if (IsValid) AddAminoAcid(nTerm.Residue);
        }

        private Tuple<double, string> GetFragmentScoreAndModifications(
            int seqIndex,
            int modIndex,
            IScorer scorer,
            IReadOnlyList<double?[]> nodeScore,
            IReadOnlyList<Tuple<double, string>[]> maxScoreAndMods)
        {
            var scoreAndMods = maxScoreAndMods[seqIndex][modIndex];
            if (scoreAndMods != null) return scoreAndMods;

            var node = _graph[seqIndex][modIndex];
            var curNodeScore = nodeScore[seqIndex][modIndex] ??
                (nodeScore[seqIndex][modIndex] = scorer.GetFragmentScore(GetComplementaryComposition(seqIndex, modIndex), GetComposition(seqIndex, modIndex)));

            var bestPrevNodeIndex = -1;
            var bestPrevNodeScore = double.NegativeInfinity;
            string bestPrevSequence = null;
            foreach (var prevNodeIndex in node.GetPrevNodeIndices())
            {
                var prevNodeScoreAndSequence = GetFragmentScoreAndModifications(seqIndex - 1, prevNodeIndex, scorer,
                    nodeScore, maxScoreAndMods);
                var prevNodeScore = prevNodeScoreAndSequence.Item1;
                if (prevNodeScore > bestPrevNodeScore)
                {
                    bestPrevNodeIndex = prevNodeIndex;
                    bestPrevNodeScore = prevNodeScore;
                    bestPrevSequence = prevNodeScoreAndSequence.Item2;
                }
            }

            if (bestPrevNodeIndex < 0)  // source
            {
                return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, string>((double)curNodeScore, "");
            }

            var modPos = _index - seqIndex;
            var aminoAcid = _aminoAcidSequence[seqIndex];
            if (aminoAcid is ModifiedAminoAcid modAa)
            {
                var modificationName = modAa.Modification.Name;
                if (string.IsNullOrEmpty(bestPrevSequence))
                {
                    bestPrevSequence = modificationName + " " + modPos;
                }
                else
                {
                    bestPrevSequence = modificationName + " " + modPos + "," + bestPrevSequence;
                }
            }

            var prevModCombIndex = _graph[seqIndex - 1][bestPrevNodeIndex].ModificationCombinationIndex;
            var curModCombIndex = node.ModificationCombinationIndex;
            if (prevModCombIndex != curModCombIndex) // modified
            {
                var modificationName = ModificationParams.GetModificationIndexBetween(prevModCombIndex, curModCombIndex).Name;
                string newModSequence;
                if (string.IsNullOrEmpty(bestPrevSequence))
                {
                    newModSequence = modificationName + " " + modPos;
                }
                else
                {
                    newModSequence = modificationName + " " + modPos + ",";
                }
                return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, string>
                    ((double)curNodeScore + bestPrevNodeScore, newModSequence+bestPrevSequence);
            }

            return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, string>
                ((double)curNodeScore + bestPrevNodeScore, bestPrevSequence);
        }

        /// <summary>
        /// Get the composition at the provided sequence index and mod index
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <returns></returns>
        protected Composition.Composition GetComposition(int seqIndex, int modIndex)
        {
            if (_nodeComposition[seqIndex][modIndex] == null)
            {
                var node = _graph[seqIndex][modIndex];
                _nodeComposition[seqIndex][modIndex] = _suffixComposition[seqIndex] +
                                  _modificationParams.GetModificationCombination(node.ModificationCombinationIndex).Composition;
            }
            return _nodeComposition[seqIndex][modIndex];
        }

        /// <summary>
        /// Get the complementary composition at the provided sequence index and mod index
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <returns></returns>
        protected Composition.Composition GetComplementaryComposition(int seqIndex, int modIndex)
        {
            if (_compNodeComposition[seqIndex][modIndex] == null)
            {
                var nodeComposition = GetComposition(seqIndex, modIndex);
                _compNodeComposition[seqIndex][modIndex] = _sinkSequenceComposition - nodeComposition;
            }
            return _compNodeComposition[seqIndex][modIndex];
        }

        /// <summary>
        /// Get the fragment score at the provided sequence index and mod index, using the provided scorer
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <param name="scorer"></param>
        /// <param name="nodeScore"></param>
        /// <param name="maxScore"></param>
        /// <returns></returns>
        protected double GetFragmentScore(int seqIndex, int modIndex, IScorer scorer, double?[][] nodeScore, double?[][] maxScore)
        {
            var score = maxScore[seqIndex][modIndex];
            if (score != null) return (double)score;

            var node = _graph[seqIndex][modIndex];
            /*double? curNodeScore;

            if (nodeScore[seqIndex][modIndex] != null)
            {
                curNodeScore = nodeScore[seqIndex][modIndex];
            }
            else
            {
                var prefixFragmentComposition = GetComplementaryComposition(seqIndex, modIndex);
                var suffixFragmentComposition = GetComposition(seqIndex, modIndex);

                if (scorer == null)
                {
                    ConsoleMsgUtils.ShowWarning("Null scorer in GetFragmentScore");
                    curNodeScore = 0;
                }
                else
                {
                    nodeScore[seqIndex][modIndex] = scorer.GetFragmentScore(prefixFragmentComposition,
                                                                            suffixFragmentComposition);
                    curNodeScore = nodeScore[seqIndex][modIndex];
                }
            }*/
            var curNodeScore = nodeScore[seqIndex][modIndex] ??
                (nodeScore[seqIndex][modIndex] = scorer.GetFragmentScore(GetComplementaryComposition(seqIndex, modIndex), GetComposition(seqIndex, modIndex)));

            var prevNodeScore = 0.0;
            if (node.GetPrevNodeIndices().Any() && seqIndex > 0)
            {
                prevNodeScore = node.GetPrevNodeIndices().Max(previousNode => GetFragmentScore(seqIndex - 1, previousNode, scorer, nodeScore, maxScore));
            }
            maxScore[seqIndex][modIndex] = curNodeScore + prevNodeScore;
            // ReSharper disable PossibleInvalidOperationException
            return (double)maxScore[seqIndex][modIndex];
            // ReSharper restore PossibleInvalidOperationException
        }

        /// <summary>
        /// Max sequence index
        /// </summary>
        protected readonly int _maxSeqIndex;

        private readonly ModificationParams _modificationParams;

        private int _index;

        /// <summary>
        /// The sequence graph itself
        /// </summary>
        protected readonly Node[][] _graph;

        /// <summary>
        /// The sequence being scored
        /// </summary>
        protected readonly string _sequence;

        private readonly AminoAcid _nTerm;
        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition.Composition[][] _nodeComposition;
        private readonly Composition.Composition[] _suffixComposition;

        private void SetNTerminalAminoAcid(AminoAcid nTerm)
        {
            _aminoAcidSequence[_aminoAcidSequence.Length - 1 - NumNTermCleavages] = nTerm;
        }

        private bool AddAminoAcid(char residue)
        {
            return PutAminoAcid(_index, residue);
        }

        /// <summary>
        /// Add an amino acid residue to this generator.
        /// </summary>
        /// <param name="index">index to add the amino acid. 0 is C-term. 1 is the C-term amino acid.</param>
        /// <param name="residue">amino acid residue to add.</param>
        /// <returns>true if residue is a valid amino acid; false otherwise.</returns>
        private bool PutAminoAcid(int index, char residue)
        {
            _index = index + 1;

            SequenceLocation? location = null;

            if(_index == 1) // C-term residue
            {
                if (residue == AminoAcid.PeptideCTerm.Residue) location = SequenceLocation.PeptideCTerm;
                else if (residue == AminoAcid.ProteinCTerm.Residue) location = SequenceLocation.ProteinCTerm;
            }
            else if (_index == _aminoAcidSequence.Length - 1 - NumNTermCleavages)   // N-term residue
            {
                if (residue == AminoAcid.PeptideNTerm.Residue) location = SequenceLocation.PeptideNTerm;
                else if (residue == AminoAcid.ProteinNTerm.Residue) location = SequenceLocation.ProteinNTerm;
            }
            else if(_index == 2) // Amino acid at the C-term
            {
                if (_aminoAcidSequence[1] == AminoAcid.PeptideCTerm) location = SequenceLocation.PeptideCTerm;
                else if (_aminoAcidSequence[1] == AminoAcid.ProteinCTerm) location = SequenceLocation.ProteinCTerm;
            }
            else if (_index == _aminoAcidSequence.Length - 2 - NumNTermCleavages) // Amino acid at the N-term
            {
                if (_aminoAcidSequence[_aminoAcidSequence.Length - 1] == AminoAcid.PeptideNTerm) location = SequenceLocation.PeptideNTerm;
                else if (_aminoAcidSequence[_aminoAcidSequence.Length - 1] == AminoAcid.ProteinNTerm) location = SequenceLocation.ProteinNTerm;
            }
            else
            {
                location = SequenceLocation.Everywhere;
            }

            if (location == null) return false;

            var loc = (SequenceLocation) location;
            var aminoAcid = AminoAcidSet.GetAminoAcid(residue, loc);
            if (aminoAcid == null) // residue is not valid
            {
                return false;
            }

            _aminoAcidSequence[_index] = aminoAcid;
            _suffixComposition[_index] = _suffixComposition[_index - 1] + aminoAcid.Composition;

            var modIndices = AminoAcidSet.GetModificationIndices(residue, loc);
            if (!modIndices.Any())  // No modification
            {
                _graph[_index] = new Node[_graph[_index - 1].Length];
                for (var i = 0; i < _graph[_index - 1].Length; i++)
                {
                    _graph[_index][i] = new Node(_graph[_index - 1][i].ModificationCombinationIndex, i);
                }
            }
            else
            {
                var modCombIndexToNodeMap = new Dictionary<int, Node>();
                for (var i = 0; i < _graph[_index - 1].Length; i++)
                {
                    var prevNodeIndex = i;
                    var prevNode = _graph[_index - 1][i];
                    var prevModCombIndex = prevNode.ModificationCombinationIndex;

                    // unmodified edge
                    if (modCombIndexToNodeMap.TryGetValue(prevModCombIndex, out var unmodifiedEdgeNode))
                    {
                        unmodifiedEdgeNode.AddPrevNodeIndex(prevNodeIndex);
                    }
                    else
                    {
                        modCombIndexToNodeMap.Add(prevModCombIndex, new Node(prevModCombIndex, prevNodeIndex));
                    }

                    // modified edge
                    foreach (var modIndex in modIndices)
                    {
                        var modCombIndex = ModificationParams.GetModificationCombinationIndex(
                                                    prevNode.ModificationCombinationIndex, modIndex);
                        if (modCombIndex < 0)   // too many modifications
                            continue;
                        if (modCombIndexToNodeMap.TryGetValue(modCombIndex, out var modifiedEdgeNode))
                        {
                            modifiedEdgeNode.AddPrevNodeIndex(prevNodeIndex);
                        }
                        else
                        {
                            modCombIndexToNodeMap.Add(modCombIndex, new Node(modCombIndex, prevNodeIndex));
                        }
                    }
                    _graph[_index] = modCombIndexToNodeMap.Values.ToArray();
                }
            }

            return true;
        }
    }

    /// <summary>
    /// A node in the Sequence graph
    /// </summary>
    public class Node
    {
        internal Node(int modificationCombinationIndex)
        {
            ModificationCombinationIndex = modificationCombinationIndex;
            _prevNodeIndices = new int[2];
            _count = 0;
        }

        internal Node(int modificationCombinationIndex, int prevNodeIndex)
            : this(modificationCombinationIndex)
        {
            AddPrevNodeIndex(prevNodeIndex);
        }

        /// <summary>
        /// The index of the modification combination
        /// </summary>
        public int ModificationCombinationIndex { get; }

        private int[] _prevNodeIndices;
        private int _count;

        internal bool AddPrevNodeIndex(int prevNodeIndex)
        {
            for (var i = 0; i < _count; i++)
            {
                if (_prevNodeIndices[i] == prevNodeIndex) return false;
            }
            if(_count >= _prevNodeIndices.Length)
            {
                Array.Resize(ref _prevNodeIndices, _prevNodeIndices.Length * 2);
            }
            _prevNodeIndices[_count++] = prevNodeIndex;

            return true;
        }

        /// <summary>
        /// Get the indices of the previous nodes
        /// </summary>
        /// <returns></returns>
        public IEnumerable<int> GetPrevNodeIndices()
        {
            for (var i = 0; i < _count; i++)
            {
                yield return _prevNodeIndices[i];
            }
        }
    }
}
