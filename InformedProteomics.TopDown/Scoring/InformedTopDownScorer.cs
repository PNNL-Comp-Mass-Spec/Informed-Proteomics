using System;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.Interfaces;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    public class InformedTopDownScorer
    {
        public InformedTopDownScorer(LcMsRun run, AminoAcidSet aaSet, int minProductCharge, int maxProductCharge, Tolerance tolerance, double ms2CorrThreshold = 0.7, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            Run = run;
            AminoAcidSet = aaSet;
            MinProductCharge = minProductCharge;
            MaxProductCharge = maxProductCharge;
            Tolerance = tolerance;
            Ms2CorrThreshold = ms2CorrThreshold;
            ActivationMethod = activationMethod;
        }

        public LcMsRun Run { get; }
        public AminoAcidSet AminoAcidSet { get; }
        public int MinProductCharge { get; }
        public int MaxProductCharge { get; }
        public Tolerance Tolerance { get; }
        public double Ms2CorrThreshold { get; }
        public ActivationMethod ActivationMethod { get; }

        public IcScores GetScores(Sequence sequence, int parentIonCharge, int ms2ScanNum)
        {
            GetCompositeScores(sequence, parentIonCharge, ms2ScanNum, out var score, out var nMatchedFragments);
            return new IcScores(nMatchedFragments, score, sequence.GetModificationString());
        }

        public IcScores GetScores(AminoAcid nTerm, string seqStr, AminoAcid cTerm, Composition composition, int charge, int ms2ScanNum)
        {
            if (!(Run.GetSpectrum(ms2ScanNum) is ProductSpectrum spec))
                return null;

            return GetScores(spec, seqStr, composition, charge, ms2ScanNum);
        }

        public IcScores GetIcScores(IInformedScorer informedScorer, IScorer scorer, string seqStr, Composition composition)
        {
            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, seqStr, AminoAcid.ProteinCTerm);
            if (seqGraph == null)
                return null;

            var bestScore = double.NegativeInfinity;
            Tuple<double, string> bestScoreAndModifications = null;
            var proteinCompositions = seqGraph.GetSequenceCompositions();

            for (var modIndex = 0; modIndex < proteinCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);
                var proteinCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();

                if (!proteinCompositionWithH2O.Equals(composition)) continue;

                var curScoreAndModifications = seqGraph.GetFragmentScoreAndModifications(scorer);
                var curScore = curScoreAndModifications.Item1;

                if (!(curScore > bestScore)) continue;

                bestScoreAndModifications = curScoreAndModifications;
                bestScore = curScore;
            }

            if (bestScoreAndModifications == null) return null;

            var modifications = bestScoreAndModifications.Item2;
            var sequence = Sequence.CreateSequence(seqStr, modifications, AminoAcidSet);
            var numMatchedFragments = informedScorer.GetNumMatchedFragments(sequence);
            var score = informedScorer.GetUserVisibleScore(sequence);

            return new IcScores(numMatchedFragments, score, modifications);
        }

        public IcScores GetScores(ProductSpectrum spec, string seqStr, Composition composition, int charge, int ms2ScanNum)
        {
            if (spec == null) return null;
            var scorer = new CompositeScorer(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, charge));
            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, seqStr, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return null;

            var bestScore = double.NegativeInfinity;
            Tuple<double, string> bestScoreAndModifications = null;
            var proteinCompositions = seqGraph.GetSequenceCompositions();

            for (var modIndex = 0; modIndex < proteinCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);
                var proteinCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();

                if (!proteinCompositionWithH2O.Equals(composition)) continue;

                var curScoreAndModifications = seqGraph.GetFragmentScoreAndModifications(scorer);
                var curScore = curScoreAndModifications.Item1;

                if (!(curScore > bestScore)) continue;

                bestScoreAndModifications = curScoreAndModifications;
                bestScore = curScore;
            }

            if (bestScoreAndModifications == null) return null;

            var modifications = bestScoreAndModifications.Item2;
            var seqObj = Sequence.CreateSequence(seqStr, modifications, AminoAcidSet);

            GetCompositeScores(seqObj, charge, ms2ScanNum, out var score, out var nMatchedFragments);
            return new IcScores(nMatchedFragments, score, modifications);
        }

        public void GetCompositeScores(Sequence sequence, int parentIonCharge, int ms2ScanNum, out double score, out int nMatchedFragments)
        {
            score = 0d;
            nMatchedFragments = 0;

            if (!(Run.GetSpectrum(ms2ScanNum) is ProductSpectrum spec))
                return;

            var preFixIonCheck = new bool[sequence.Count + 1];
            var sufFixIonCheck = new bool[sequence.Count + 1];

            var scorer = new CompositeScorer(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, parentIonCharge + 2));
            var cleavages = sequence.GetInternalCleavages();
            var cleavageIndex = 0;

            foreach (var c in cleavages)
            {
                score += scorer.GetFragmentScore(c.PrefixComposition, c.SuffixComposition, out var prefixHit, out var suffixHit);

                nMatchedFragments += (prefixHit) ? 1 : 0;
                nMatchedFragments += (suffixHit) ? 1 : 0;

                if (prefixHit) preFixIonCheck[cleavageIndex] = true;
                if (suffixHit) sufFixIonCheck[cleavageIndex] = true;

                cleavageIndex++;
            }

            var preContCount = 0;
            var sufContCount = 0;
            for (var i = 0; i < preFixIonCheck.Length - 1; i++)
            {
                if (preFixIonCheck[i] && preFixIonCheck[i + 1]) preContCount++;
                if (sufFixIonCheck[i] && sufFixIonCheck[i + 1]) sufContCount++;
            }
            score += preContCount * CompositeScorer.ScoreParam.Prefix.ConsecutiveMatch;
            score += sufContCount * CompositeScorer.ScoreParam.Suffix.ConsecutiveMatch;
        }

        public void GetCompositeScores(Sequence sequence, CompositeScorerBasedOnDeconvolutedSpectrum scorer, out double score)
        {
            score = 0d;

            //var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            //if (spec == null) return;

            //var scorer = new CompositeScorer(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, parentIonCharge + 2), activationMethod: ActivationMethod);
            var cleavages = sequence.GetInternalCleavages();

            foreach (var c in cleavages)
            {
                score += scorer.GetFragmentScore(c.PrefixComposition, c.SuffixComposition);
            }
        }
    }
}
