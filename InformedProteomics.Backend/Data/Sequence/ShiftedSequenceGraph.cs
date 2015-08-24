using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class ShiftedSequenceGraph
    {
        public ShiftedSequenceGraph(AminoAcidSet aminoAcidSet, double shiftedMass, bool isForward, int maxSequenceLength, double maxSequenceMass = 50000.0)
        {
            _aminoAcidSet = aminoAcidSet;
            _modificationParams = aminoAcidSet.GetModificationParams();

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
                _compNodeComposition[i] = new Composition.Composition[_modificationParams.NumModificationCombinations];
                _nodeComposition[i] = new Composition.Composition[_modificationParams.NumModificationCombinations];
            }

            IsValid = true;
        }

        public bool AddAminoAcid(char residue, SequenceLocation location = SequenceLocation.Everywhere)
        {
            return PutAminoAcid(_index, residue, location);
        }

        public AminoAcidSet AminoAcidSet
        {
            get { return _aminoAcidSet; }
        }
        public ModificationParams ModificationParams
        {
            get { return _modificationParams; }
        }

        public double ShiftMass { get; private set; }

        public bool IsValid { get; private set; }


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

        public void SetSink(int modIndex)
        {
            _sinkModIndex = modIndex;
            _sinkSequenceComposition = GetComposition(_index, modIndex);
            _sinkSequenceCompositionWithH2O = _sinkSequenceComposition + Composition.Composition.H2O;
            //_compNodeComposition = new Composition.Composition[_maxSeqIndex, _modificationParams.NumModificationCombinations];
            foreach (var t in _compNodeComposition) Array.Clear(t, 0, t.Length);
        }

        public Composition.Composition GetSinkSequenceCompositionWithH2O()
        {
            return _sinkSequenceCompositionWithH2O;
        }
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

        private Tuple<double, LinkedList<ModificationInstance>> GetFragmentScoreAndModifications(int seqIndex, int modIndex,
            IScorer scorer, double?[][] nodeScore, Tuple<double, LinkedList<ModificationInstance>>[][] maxScoreAndMods)
        {
            var scoreAndMods = maxScoreAndMods[seqIndex][modIndex];
            if (scoreAndMods != null) return scoreAndMods;

            var node = _graph[seqIndex][modIndex];
            var curNodeScore = nodeScore[seqIndex][modIndex] ??
                               (nodeScore[seqIndex][modIndex] = _isForward
                                   ? scorer.GetFragmentScore(GetComposition(seqIndex, modIndex),
                                       GetComplementaryComposition(seqIndex, modIndex))
                                   : scorer.GetFragmentScore(GetComplementaryComposition(seqIndex, modIndex),
                                       GetComposition(seqIndex, modIndex))
                                   );

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
            if (modPos <= 1) --modPos;
            var aminoAcid = _aminoAcidSequence[seqIndex];
            var modAa = aminoAcid as ModifiedAminoAcid;

            LinkedList<ModificationInstance> newMods = null;
            if (modAa != null)
            {
                newMods = bestPrevMods == null ? new LinkedList<ModificationInstance>() : new LinkedList<ModificationInstance>(bestPrevMods);
                var modIns = new ModificationInstance(modAa.Modification, modPos);
                if (_isForward) newMods.AddLast(modIns);
                else newMods.AddFirst(modIns);
            }

            var prevModCombIndex = _graph[seqIndex - 1][bestPrevNodeIndex].ModificationCombinationIndex;
            var curModCombIndex = node.ModificationCombinationIndex;
            if (prevModCombIndex != curModCombIndex) // modified
            {
                if (newMods == null)
                {
                    newMods = bestPrevMods == null ? new LinkedList<ModificationInstance>() : new LinkedList<ModificationInstance>(bestPrevMods);
                }
                var modification = ModificationParams.GetModificationIndexBetween(prevModCombIndex, curModCombIndex);
                var modIns = new ModificationInstance(modification, modPos);
                if (_isForward) newMods.AddLast(modIns);
                else newMods.AddFirst(modIns);
            }

            return maxScoreAndMods[seqIndex][modIndex] = new Tuple<double, LinkedList<ModificationInstance>>
                ((double)curNodeScore + bestPrevNodeScore, newMods ?? bestPrevMods);
        }

        protected Composition.Composition GetComposition(int seqIndex, int modIndex)
        {
            if (_nodeComposition[seqIndex][modIndex] == null)
            {
                var node = _graph[seqIndex][modIndex];
                _nodeComposition[seqIndex][modIndex] = _fragmentComposition[seqIndex] +
                                  _modificationParams.GetModificationCombination(node.ModificationCombinationIndex)
                                                     .Composition;
            }
            return _nodeComposition[seqIndex][modIndex];
        }

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

        protected double GetFragmentScore(int seqIndex, int modIndex, IScorer scorer, double?[][] nodeScore, double?[][] maxScore)
        {
            var score = maxScore[seqIndex][modIndex];
            if (score != null) return (double)score;

            var node = _graph[seqIndex][modIndex];
            var curNodeScore = nodeScore[seqIndex][modIndex] ??
                (nodeScore[seqIndex][modIndex] = 
                    _isForward ? scorer.GetFragmentScore(GetComposition(seqIndex, modIndex), GetComplementaryComposition(seqIndex, modIndex))
                    : scorer.GetFragmentScore(GetComplementaryComposition(seqIndex, modIndex), GetComposition(seqIndex, modIndex))
                    );

            var prevNodeScore = 0.0;
            if (node.GetPrevNodeIndices().Any())
            {
                prevNodeScore =
                    node.GetPrevNodeIndices()
                        .Max(prevNodeIndex => GetFragmentScore(seqIndex - 1, prevNodeIndex, scorer, nodeScore, maxScore));
            }
            maxScore[seqIndex][modIndex] = curNodeScore + prevNodeScore;
            // ReSharper disable PossibleInvalidOperationException
            return (double)maxScore[seqIndex][modIndex];
            // ReSharper restore PossibleInvalidOperationException
        }

        private Composition.Composition _sinkSequenceComposition;
        private Composition.Composition _sinkSequenceCompositionWithH2O;
        private int _sinkModIndex;
        //private readonly Composition.Composition[,] _compNodeComposition;
        private readonly Composition.Composition[][] _compNodeComposition;

        private readonly int _maxSeqIndex;
        private readonly AminoAcidSet _aminoAcidSet;
        private readonly ModificationParams _modificationParams;
        private readonly bool _isForward;
        private readonly double _maxSequenceMass;

        private int _index;
        private readonly Node[][] _graph;

        private readonly AminoAcid[] _aminoAcidSequence;
        private readonly Composition.Composition[][] _nodeComposition;

        private readonly Composition.Composition[] _fragmentComposition;

        /// <summary>
        /// Add an amino acid residue to this generator.
        /// </summary>
        /// <param name="index">index to add the amino acid. 0 is C-term. 1 is the C-term amino acid.</param>
        /// <param name="residue">amino acid residue to add.</param>
        /// <param name="loc">location of the residue</param>
        /// <returns>true if residue is a valid amino acid; false otherwise.</returns>
        private bool PutAminoAcid(int index, char residue, SequenceLocation loc)
        {
            _index = index + 1;

            var aminoAcid = AminoAcidSet.GetAminoAcid(residue, loc);
            if (aminoAcid == null) // residue is not valid
            {
                return false;
            }

            var fragmentComposition = _fragmentComposition[_index - 1] + aminoAcid.Composition;
            if (fragmentComposition.Mass > _maxSequenceMass) return false;

            _aminoAcidSequence[_index] = aminoAcid;
            _fragmentComposition[_index] = fragmentComposition;

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
                    Node newNode;
                    // unmodified edge
                    if (modCombIndexToNodeMap.TryGetValue(prevModCombIndex, out newNode))
                    {
                        newNode.AddPrevNodeIndex(prevNodeIndex);
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
                        if (modCombIndexToNodeMap.TryGetValue(modCombIndex, out newNode))
                        {
                            newNode.AddPrevNodeIndex(prevNodeIndex);
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
