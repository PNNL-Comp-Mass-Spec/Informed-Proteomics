using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// A shifted mass sequence graph
    /// </summary>
    // ReSharper disable once UnusedMember.Global
    public class ShiftedSequenceGraph
    {
        // Ignore Spelling: proteoforms

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="aminoAcidSet"></param>
        /// <param name="shiftedMass"></param>
        /// <param name="isForward"></param>
        /// <param name="maxSequenceLength"></param>
        /// <param name="maxSequenceMass"></param>
        public ShiftedSequenceGraph(AminoAcidSet aminoAcidSet, double shiftedMass, bool isForward, int maxSequenceLength, double maxSequenceMass = 50000.0)
        {
            AminoAcidSet = aminoAcidSet;
            ModificationParams = aminoAcidSet.GetModificationParams();

            _isForward = isForward;

            _index = 0;
            _maxSeqIndex = maxSequenceLength + 2;   // shift + Term + length
            _maxSequenceMass = maxSequenceMass;

            _aminoAcidSequence = new AminoAcid[_maxSeqIndex];
            var shiftAa = new AminoAcid('\0', "Shift", new CompositionWithDeltaMass(shiftedMass));
            _aminoAcidSequence[0] = shiftAa;

            ShiftMass = shiftedMass;

            _fragmentComposition = new Composition.Composition[_maxSeqIndex];
            _fragmentComposition[0] = shiftAa.Composition;

            _graph = new Node[_maxSeqIndex][];
            _graph[0] = new[] { new Node(0) };

            _nodeComposition = new Composition.Composition[_maxSeqIndex][];
            _compNodeComposition = new Composition.Composition[_maxSeqIndex][];
            for (var i = 0; i < _maxSeqIndex; i++)
            {
                _compNodeComposition[i] = new Composition.Composition[ModificationParams.NumModificationCombinations];
                _nodeComposition[i] = new Composition.Composition[ModificationParams.NumModificationCombinations];
            }

            IsValid = true;
        }

        /// <summary>
        /// Add an amino acid to this shifted sequence graph
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="location"></param>
        /// <returns>True if residue is a valid amino acid, otherwise false</returns>
        public bool AddAminoAcid(char residue, SequenceLocation location = SequenceLocation.Everywhere)
        {
            return PutAminoAcid(_index, residue, location);
        }

        /// <summary>
        /// The amino acid set used in this shifted sequence graph
        /// </summary>
        public AminoAcidSet AminoAcidSet { get; }

        /// <summary>
        /// Modifications used in this shifted sequence graph
        /// </summary>
        public ModificationParams ModificationParams { get; }

        /// <summary>
        /// The mass by which this sequence graph is shifted
        /// </summary>
        public double ShiftMass { get; }

        /// <summary>
        /// If this shifted sequence graph is valid
        /// </summary>
        public bool IsValid { get; }

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
        /// Set the sink to the provided modification index
        /// </summary>
        /// <param name="modIndex"></param>
        public void SetSink(int modIndex)
        {
            _sinkModIndex = modIndex;
            _sinkSequenceComposition = GetComposition(_index, modIndex);
            _sinkSequenceCompositionWithH2O = _sinkSequenceComposition + Composition.Composition.H2O;
            //_compNodeComposition = new Composition.Composition[_maxSeqIndex, _modificationParams.NumModificationCombinations];
            foreach (var t in _compNodeComposition)
            {
                Array.Clear(t, 0, t.Length);
            }
        }

        /// <summary>
        /// Get the sink sequence composition with H2O added
        /// </summary>
        /// <returns>Composition object</returns>
        public Composition.Composition GetSinkSequenceCompositionWithH2O()
        {
            return _sinkSequenceCompositionWithH2O;
        }

        /// <summary>
        /// Get the fragment score using the provided scorer
        /// </summary>
        /// <param name="scorer"></param>
        /// <returns>Fragment score</returns>
        public double GetFragmentScore(IScorer scorer)
        {
            var nodeScore = new double?[_maxSeqIndex][];
            var maxScore = new double?[_maxSeqIndex][];
            for (var si = 0; si <= _index; si++)
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
        /// Get the score and modifications
        /// </summary>
        /// <param name="scorer"></param>
        /// <returns>Tuple of fragment score and list of modifications</returns>
        public Tuple<double, LinkedList<ModificationInstance>> GetScoreAndModifications(IScorer scorer)
        {
            var nodeScore = new double?[_maxSeqIndex][];
            var maxScore = new Tuple<double, LinkedList<ModificationInstance>>[_maxSeqIndex][];
            for (var si = 0; si <= _index; si++)
            {
                nodeScore[si] = new double?[_graph[si].Length];
                maxScore[si] = new Tuple<double, LinkedList<ModificationInstance>>[_graph[si].Length];
            }
            maxScore[0][0] = new Tuple<double, LinkedList<ModificationInstance>>(0.0, new LinkedList<ModificationInstance>());
            nodeScore[_index][_sinkModIndex] = 0.0;

            var fragmentScore = GetFragmentScoreAndModifications(_index, _sinkModIndex, scorer, nodeScore,
                maxScore);

            return new Tuple<double, LinkedList<ModificationInstance>>(fragmentScore.Item1, fragmentScore.Item2);
        }

        /// <summary>
        /// Get the fragment score and modifications
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <param name="scorer"></param>
        /// <param name="nodeScore"></param>
        /// <param name="maxScoreAndMods"></param>
        /// <returns>Tuple of fragment score and list of modifications</returns>
        private Tuple<double, LinkedList<ModificationInstance>> GetFragmentScoreAndModifications(
            int seqIndex,
            int modIndex,
            IScorer scorer,
            IReadOnlyList<double?[]> nodeScore,
            IReadOnlyList<Tuple<double, LinkedList<ModificationInstance>>[]> maxScoreAndMods)
        {
            var scoreAndMods = maxScoreAndMods[seqIndex][modIndex];
            if (scoreAndMods != null)
            {
                return scoreAndMods;
            }

            var node = _graph[seqIndex][modIndex];

            // The following will use nodeScore[seqIndex][modIndex] if it is not null
            // Otherwise, it will update nodeScore[seqIndex][modIndex] using scorer.GetFragmentScore(), then store in curNodeScore
            var curNodeScore =
                    nodeScore[seqIndex][modIndex] ??=
                        _isForward
                            ? scorer.GetFragmentScore(
                                GetComposition(seqIndex, modIndex),
                                GetComplementaryComposition(seqIndex, modIndex))
                            : scorer.GetFragmentScore(
                                GetComplementaryComposition(seqIndex, modIndex),
                                GetComposition(seqIndex, modIndex))
;

            var bestPrevNodeIndex = -1;
            var bestPrevNodeScore = double.NegativeInfinity;
            LinkedList<ModificationInstance> bestPrevMods = null;
            foreach (var prevNodeIndex in node.GetPrevNodeIndices())
            {
                var prevNodeScoreAndSequence = GetFragmentScoreAndModifications(seqIndex - 1, prevNodeIndex, scorer,
                    nodeScore, maxScoreAndMods);
                var prevNodeScore = prevNodeScoreAndSequence.Item1;
                if (prevNodeScore > bestPrevNodeScore)
                {
                    bestPrevNodeIndex = prevNodeIndex;
                    bestPrevNodeScore = prevNodeScore;
                    bestPrevMods = prevNodeScoreAndSequence.Item2;
                }
            }

            if (bestPrevNodeIndex < 0)  // source
            {
                return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, LinkedList<ModificationInstance>>((double)curNodeScore, new LinkedList<ModificationInstance>());
            }

            var modPos = _isForward ? seqIndex : _index - seqIndex;
            if (modPos <= 1)
            {
                --modPos;
            }

            var aminoAcid = _aminoAcidSequence[seqIndex];

            LinkedList<ModificationInstance> newMods = null;
            if (aminoAcid is ModifiedAminoAcid modAa)
            {
                newMods = bestPrevMods == null ? new LinkedList<ModificationInstance>() : new LinkedList<ModificationInstance>(bestPrevMods);
                var modIns = new ModificationInstance(modAa.Modification, modPos);
                if (_isForward)
                {
                    newMods.AddLast(modIns);
                }
                else
                {
                    newMods.AddFirst(modIns);
                }
            }

            var prevModCombIndex = _graph[seqIndex - 1][bestPrevNodeIndex].ModificationCombinationIndex;
            var curModCombIndex = node.ModificationCombinationIndex;
            if (prevModCombIndex != curModCombIndex) // modified
            {
                newMods ??= bestPrevMods == null ? new LinkedList<ModificationInstance>() : new LinkedList<ModificationInstance>(bestPrevMods);

                var modification = ModificationParams.GetModificationIndexBetween(prevModCombIndex, curModCombIndex);
                var modIns = new ModificationInstance(modification, modPos);
                if (_isForward)
                {
                    newMods.AddLast(modIns);
                }
                else
                {
                    newMods.AddFirst(modIns);
                }
            }

            return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, LinkedList<ModificationInstance>>
                ((double)curNodeScore + bestPrevNodeScore, newMods ?? bestPrevMods);
        }

        /// <summary>
        /// Get the composition at the provided sequence index and mod index
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <returns>Composition object</returns>
        protected Composition.Composition GetComposition(int seqIndex, int modIndex)
        {
            if (_nodeComposition[seqIndex][modIndex] == null)
            {
                var node = _graph[seqIndex][modIndex];
                _nodeComposition[seqIndex][modIndex] = _fragmentComposition[seqIndex] +
                                  ModificationParams.GetModificationCombination(node.ModificationCombinationIndex)
                                                     .Composition;
            }
            return _nodeComposition[seqIndex][modIndex];
        }

        /// <summary>
        /// Get the complementary composition at the provided sequence index and mod index
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <returns>Composition object</returns>
        protected Composition.Composition GetComplementaryComposition(int seqIndex, int modIndex)
        {
            if (_compNodeComposition[seqIndex][modIndex] == null)
            {
                var nodeComposition = GetComposition(seqIndex, modIndex);
                _compNodeComposition[seqIndex][modIndex] = _sinkSequenceComposition - nodeComposition;
            }
            return _compNodeComposition[seqIndex][modIndex];
            /*
            if (_compNodeComposition[seqIndex, modIndex] == null)
            {
                var nodeComposition = GetComposition(seqIndex, modIndex);
                _compNodeComposition[seqIndex, modIndex] = _sinkSequenceComposition - nodeComposition;
            }
            return _compNodeComposition[seqIndex, modIndex];*/
        }

        /// <summary>
        /// Get the fragment score at the provided sequence index and mod index, using the provided scorer
        /// </summary>
        /// <param name="seqIndex"></param>
        /// <param name="modIndex"></param>
        /// <param name="scorer"></param>
        /// <param name="nodeScore"></param>
        /// <param name="maxScore"></param>
        /// <returns>Fragment score</returns>
        protected double GetFragmentScore(int seqIndex, int modIndex, IScorer scorer, double?[][] nodeScore, double?[][] maxScore)
        {
            var score = maxScore[seqIndex][modIndex];
            if (score != null)
            {
                return (double)score;
            }

            var node = _graph[seqIndex][modIndex];

            // The following will use nodeScore[seqIndex][modIndex] if it is not null
            // Otherwise, it will update nodeScore[seqIndex][modIndex] using scorer.GetFragmentScore(), then store in curNodeScore
            var curNodeScore =
                    nodeScore[seqIndex][modIndex] ??=
                        _isForward
                        ? scorer.GetFragmentScore(
                            GetComposition(seqIndex, modIndex),
                            GetComplementaryComposition(seqIndex, modIndex))
                        : scorer.GetFragmentScore(
                            GetComplementaryComposition(seqIndex, modIndex),
                            GetComposition(seqIndex, modIndex));

            var prevNodeScore = 0.0;
            if (node.GetPrevNodeIndices().Any())
            {
                prevNodeScore =
                    node.GetPrevNodeIndices()
                        .Max(prevNodeIndex => GetFragmentScore(seqIndex - 1, prevNodeIndex, scorer, nodeScore, maxScore));
            }
            maxScore[seqIndex][modIndex] = curNodeScore + prevNodeScore;

            // ReSharper disable once PossibleInvalidOperationException
            return (double)maxScore[seqIndex][modIndex];
        }

        private Composition.Composition _sinkSequenceComposition;
        private Composition.Composition _sinkSequenceCompositionWithH2O;
        private int _sinkModIndex;

        private readonly Composition.Composition[][] _compNodeComposition;

        private readonly int _maxSeqIndex;
        private readonly bool _isForward;
        private readonly double _maxSequenceMass;

        private int _index;
        private readonly Node[][] _graph;

        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition.Composition[][] _nodeComposition;

        private readonly Composition.Composition[] _fragmentComposition;

        /// <summary>
        /// Add an amino acid residue to this generator
        /// </summary>
        /// <param name="index">Index to add the amino acid; 0 is C-terminus; 1 is the C-terminal amino acid</param>
        /// <param name="residue">Amino acid residue to add</param>
        /// <param name="loc">Location of the residue</param>
        /// <returns>True if residue is a valid amino acid, otherwise false</returns>
        private bool PutAminoAcid(int index, char residue, SequenceLocation loc)
        {
            _index = index + 1;

            var aminoAcid = AminoAcidSet.GetAminoAcid(residue, loc);
            if (aminoAcid == null) // residue is not valid
            {
                return false;
            }

            var fragmentComposition = _fragmentComposition[_index - 1] + aminoAcid.Composition;
            if (fragmentComposition.Mass > _maxSequenceMass)
            {
                return false;
            }

            _aminoAcidSequence[_index] = aminoAcid;
            _fragmentComposition[_index] = fragmentComposition;

            var modIndices = AminoAcidSet.GetModificationIndices(residue, loc);
            if (modIndices.Length == 0)  // No modification
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
                        {
                            continue;
                        }

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
}
