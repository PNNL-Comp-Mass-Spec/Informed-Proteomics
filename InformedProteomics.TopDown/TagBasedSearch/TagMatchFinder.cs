using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.TopDown.PostProcessing;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class TagMatchFinder
    {
        // TODO: better to compute this from mod file
        public const double MinSumModificationMasses = -200.0;
        public const double Ms1CorrThreshold = 0.7;

        public TagMatchFinder(
            ProductSpectrum spec,
            IScorer ms2Scorer,
            LcMsPeakMatrix featureFinder,
            string proteinSequence, 
            Tolerance tolerance, 
            AminoAcidSet aaSet, 
            double maxSequenceMass)
        {
            _spec = spec;
            _ms2Scorer = ms2Scorer;
            _featureFinder = featureFinder;
            _proteinSequence = proteinSequence;
            _tolerance = tolerance;
            _aaSet = aaSet;
            _maxSequenceMass = maxSequenceMass;
        }

        public IEnumerable<TagMatch> FindMatches(MatchedTag matchedTag)
        {
            if(matchedTag.NTermFlankingMass != null && matchedTag.CTermFlankingMass != null) return FindMatchesWithFeatureMass(matchedTag);
            if(matchedTag.NTermFlankingMass != null) return FindMatchesForwardAndBackward(matchedTag);
            if(matchedTag.CTermFlankingMass != null) return FindMatchesBackwardAndForward(matchedTag);
            return Enumerable.Empty<TagMatch>();
        }

        private IEnumerable<TagMatch> FindMatchesWithFeatureMass(MatchedTag matchedTag)
        {
            if (matchedTag.NTermFlankingMass == null || matchedTag.CTermFlankingMass == null) yield break;
            var featureMass = (double) matchedTag.NTermFlankingMass + matchedTag.Mass +
                              (double)matchedTag.CTermFlankingMass + Composition.H2O.Mass;
            var shiftMass = matchedTag.Mass + (double)matchedTag.NTermFlankingMass;

            var backwardGraph = new ShiftedSequenceGraph(_aaSet, shiftMass, false,
                matchedTag.StartIndex, featureMass - MinSumModificationMasses);

            foreach (var backwardMatch in GetBackwardMatches(matchedTag, backwardGraph, featureMass))
            {
                // Make a forward graph
                var nTermShiftMass = backwardMatch.Mass + matchedTag.Mass;
                var forwardGraph = new ShiftedSequenceGraph(_aaSet, nTermShiftMass, true,
                    _proteinSequence.Length - matchedTag.EndIndex, featureMass - MinSumModificationMasses);

                foreach (
                    var forwardMatch in
                        GetForwardMatches(matchedTag, forwardGraph, featureMass))
                {
                    var mass = forwardMatch.Mass + matchedTag.Mass + backwardMatch.Mass;
                    if (mass > _maxSequenceMass) continue;

                    var offset = matchedTag.EndIndex - backwardMatch.Index - 1;
                    var modStr = string.Join(",", backwardMatch.Modifications.Concat(forwardMatch.Modifications.Select(m => m.GetModificationInstanceWithOffset(offset))));

                    var modList = new List<Modification>();
                    foreach (var mod in backwardMatch.Modifications) modList.Add(mod.Modification);
                    foreach (var mod in forwardMatch.Modifications) modList.Add(mod.Modification);

                    var tagMatch = new TagMatch(
                        backwardMatch.Index, 
                        forwardMatch.Index, 
                        matchedTag.Length,
                        backwardMatch.Charge,
                        backwardMatch.Score,
                        forwardMatch.Score,
                        mass,
                        new ModificationCombination(modList), 
                        modStr);
                    yield return tagMatch;
                }
            }
        }

        private IEnumerable<TagMatch> FindMatchesForwardAndBackward(MatchedTag matchedTag)
        {
            if (matchedTag.NTermFlankingMass == null) yield break;
            var shiftMass = matchedTag.Mass + (double)matchedTag.NTermFlankingMass;

            var forwardGraph = new ShiftedSequenceGraph(_aaSet, shiftMass, true,
                _proteinSequence.Length - matchedTag.EndIndex, _maxSequenceMass);

            foreach (var forwardMatch in GetForwardMatches(matchedTag, forwardGraph))
            {
                // Make a backward graph
                var cTermShiftMass = forwardMatch.Mass + matchedTag.Mass;
                var featureMass = (double)matchedTag.NTermFlankingMass + cTermShiftMass + Composition.H2O.Mass;
                var backwardGraph = new ShiftedSequenceGraph(_aaSet, cTermShiftMass, false,
                    matchedTag.StartIndex, featureMass - MinSumModificationMasses);

                foreach (
                    var backwardMatch in
                        GetBackwardMatches(matchedTag, backwardGraph, featureMass))
                {
                    var mass = forwardMatch.Mass + matchedTag.Mass + backwardMatch.Mass;
                    if (mass > _maxSequenceMass) continue;

                    var offset = matchedTag.EndIndex - backwardMatch.Index - 1;
                    var modStr = string.Join(",", backwardMatch.Modifications.Concat(forwardMatch.Modifications.Select(m => m.GetModificationInstanceWithOffset(offset))));
                    
                    var modList = new List<Modification>();
                    foreach (var mod in backwardMatch.Modifications) modList.Add(mod.Modification);
                    foreach (var mod in forwardMatch.Modifications) modList.Add(mod.Modification);

                    var tagMatch = new TagMatch(
                        backwardMatch.Index, 
                        forwardMatch.Index,
                        matchedTag.Length,
                        forwardMatch.Charge,
                        backwardMatch.Score,
                        forwardMatch.Score,
                        mass,
                        new ModificationCombination(modList), 
                        modStr);
                    yield return tagMatch;
                }
            }
        }

        private IEnumerable<TagMatch> FindMatchesBackwardAndForward(MatchedTag matchedTag)
        {
            if (matchedTag.CTermFlankingMass == null) yield break;
            var shiftMass = matchedTag.Mass + (double)matchedTag.CTermFlankingMass;

            var backwardGraph = new ShiftedSequenceGraph(_aaSet, shiftMass, false,
                matchedTag.StartIndex, _maxSequenceMass);

            foreach (var backwardMatch in GetBackwardMatches(matchedTag, backwardGraph))
            {
                // Make a forward graph
                var nTermShiftMass = backwardMatch.Mass + matchedTag.Mass;
                var featureMass = nTermShiftMass + (double)matchedTag.CTermFlankingMass + Composition.H2O.Mass;
                var forwardGraph = new ShiftedSequenceGraph(_aaSet, nTermShiftMass, true,
                    _proteinSequence.Length - matchedTag.EndIndex, featureMass - MinSumModificationMasses);

                foreach (
                    var forwardMatch in
                        GetForwardMatches(matchedTag, forwardGraph, featureMass))
                {
                    var mass = forwardMatch.Mass + matchedTag.Mass + backwardMatch.Mass;
                    if (mass > _maxSequenceMass) continue;

                    var offset = matchedTag.EndIndex - backwardMatch.Index - 1;
                    var modStr = string.Join(",", backwardMatch.Modifications.Concat(forwardMatch.Modifications.Select(m => m.GetModificationInstanceWithOffset(offset))));
                    var modList = new List<Modification>();
                    foreach (var mod in backwardMatch.Modifications) modList.Add(mod.Modification);
                    foreach (var mod in forwardMatch.Modifications) modList.Add(mod.Modification);

                    var tagMatch = new TagMatch(
                        backwardMatch.Index, 
                        forwardMatch.Index,
                        matchedTag.Length,
                        backwardMatch.Charge,
                        backwardMatch.Score,
                        forwardMatch.Score,
                        mass,
                        new ModificationCombination(modList), 
                        modStr);
                    yield return tagMatch;
                }
            }
        }

        private readonly ProductSpectrum _spec;
        private readonly IScorer _ms2Scorer;
        
        private readonly LcMsPeakMatrix _featureFinder;
        private readonly string _proteinSequence;

        private readonly Tolerance _tolerance;
        private readonly AminoAcidSet _aaSet;

        private readonly double _maxSequenceMass;
//        private readonly int _minProductIonCharge;
//        private readonly int _maxProductIonCharge;

        private IEnumerable<FlankingMassMatch> GetForwardMatches(
            MatchedTag matchedTag,
            ShiftedSequenceGraph forwardGraph,
            double? featureMass = null
            )
        {
            for (var i = matchedTag.EndIndex; i <= _proteinSequence.Length; i++)
            {
                var residue = i < _proteinSequence.Length ? _proteinSequence[i] : AminoAcid.ProteinCTerm.Residue;
                var location = i < _proteinSequence.Length - 1
                    ? SequenceLocation.Everywhere
                    : SequenceLocation.ProteinCTerm;
                if (!forwardGraph.AddAminoAcid(residue, location)) yield break;

                if (i == _proteinSequence.Length - 1) continue;

                var forwardMatch = GetBestMatchInTheGraph(forwardGraph, _spec, featureMass);

                if (forwardMatch != null)
                {
                    forwardMatch.Index = Math.Min(i + 1, _proteinSequence.Length);
                    yield return forwardMatch;
                }
            }
        }

        private IEnumerable<FlankingMassMatch> GetBackwardMatches(
            MatchedTag matchedTag,
            ShiftedSequenceGraph backwardGraph,
            double? featureMass = null
            )
        {
            for (var j = matchedTag.StartIndex - 1; j >= -1; j--)
            {
                var residue = j >= 0 ? _proteinSequence[j] : AminoAcid.ProteinNTerm.Residue;
                var location = j > 0 ? SequenceLocation.Everywhere : SequenceLocation.ProteinNTerm;
                if(!backwardGraph.AddAminoAcid(residue, location)) yield break;

                if (j == 0) continue;
                var backwardMatch = GetBestMatchInTheGraph(backwardGraph, _spec, featureMass);
                if (backwardMatch != null)
                {
                    backwardMatch.Index = Math.Max(j, 0);
                    yield return backwardMatch;
                }
            }
        }

        private FlankingMassMatch GetBestMatchInTheGraph(ShiftedSequenceGraph seqGraph, ProductSpectrum spec, double? featureMass)
        {
            FlankingMassMatch match = null;
            var bestScore = double.NegativeInfinity;
            var protCompositions = seqGraph.GetSequenceCompositions();
            for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);
                var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                var sequenceMass = protCompositionWithH2O.Mass;
                if (featureMass != null && !_tolerance.IsWithin(sequenceMass, (double)featureMass)) continue;

                var charge =
                    (int)
                        Math.Round(sequenceMass /
                                   (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));

                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(sequenceMass, charge,
                    Averagine.GetIsotopomerEnvelope(sequenceMass).MostAbundantIsotopeIndex);

                if (!spec.IsolationWindow.Contains(mostAbundantIsotopeMz)) continue;

                //var feature = new TargetFeature(sequenceMass, charge, spec.ScanNum);
                
                if (_featureFinder != null)
                {
                    var ms1Corr = _featureFinder.GetMs1EvidenceScore(spec.ScanNum, sequenceMass, charge);
                    if (ms1Corr < Ms1CorrThreshold) continue;
                }

                var curScoreAndModifications = seqGraph.GetScoreAndModifications(_ms2Scorer);
                var curScore = curScoreAndModifications.Item1;
//                var curScore = seqGraph.GetFragmentScore(_ms2Scorer);
                if (curScore > bestScore)
                {
                    match = new FlankingMassMatch(curScore,
                        sequenceMass - Composition.H2O.Mass - seqGraph.ShiftMass, charge, curScoreAndModifications.Item2);
                    //match = new FlankingMassMatch(curScore,
                    //    sequenceMass - Composition.H2O.Mass - seqGraph.ShiftMass, charge, new ModificationInstance[0]);
                    bestScore = curScore;
                }
            }

            return match;
        }
    }

    internal class FlankingMassMatch
    {
        internal FlankingMassMatch(double score, double mass, int charge, IEnumerable<ModificationInstance> modifications)
        {
            Score = score;
            Mass = mass;
            Charge = charge;
            Modifications = modifications;
        }

        internal double Score { get; set; }
        internal double Mass { get; set; }
        internal int Charge { get; set; }
        internal IEnumerable<ModificationInstance> Modifications { get; set; }
        internal int Index { get; set; }

        public override string ToString()
        {
            return string.Format("{0}+, {1:F1} Da, Score {2:F1}", Charge, Mass, Score);
        }
    }
}
