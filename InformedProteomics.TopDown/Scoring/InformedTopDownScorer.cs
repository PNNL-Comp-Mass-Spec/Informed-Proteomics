using System;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

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

        public LcMsRun Run { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinProductCharge { get; private set; }
        public int MaxProductCharge { get; private set; }
        public Tolerance Tolerance { get; private set; }
        public double Ms2CorrThreshold { get; private set; }
        public ActivationMethod ActivationMethod { get; private set; }

        public IcScores GetScores(Sequence sequence, int parentIoncharge, int ms2ScanNum)
        {
            double score;
            int nMatchedFragments;
            GetCompositeScores(sequence, parentIoncharge, ms2ScanNum, out score, out nMatchedFragments);
            return new IcScores(nMatchedFragments, score, sequence.GetModificationString());
        }

        public IcScores GetScores(AminoAcid nTerm, string seqStr, AminoAcid cTerm, Composition composition, int charge, int ms2ScanNum)
        {
            var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return null;
            return GetScores(spec, seqStr, composition, charge, ms2ScanNum);
        }

        public IcScores GetScores(ProductSpectrum spec, string seqStr, Composition composition, int charge, int ms2ScanNum)
        {
            if (spec == null) return null;
            var scorer = new CompositeScorer(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, charge), activationMethod: ActivationMethod);
            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, seqStr, AminoAcid.ProteinCTerm);
            if (seqGraph == null) return null;

            var bestScore = double.NegativeInfinity;
            Tuple<double, string> bestScoreAndModifications = null;
            var protCompositions = seqGraph.GetSequenceCompositions();

            for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);
                var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();

                if (!protCompositionWithH2O.Equals(composition)) continue;

                var curScoreAndModifications = seqGraph.GetFragmentScoreAndModifications(scorer);
                var curScore = curScoreAndModifications.Item1;

                if (!(curScore > bestScore)) continue;

                bestScoreAndModifications = curScoreAndModifications;
                bestScore = curScore;
            }

            if (bestScoreAndModifications == null) return null;

            var modifications = bestScoreAndModifications.Item2;
            var seqObj = Sequence.CreateSequence(seqStr, modifications, AminoAcidSet);

            double score;
            int nMatchedFragments;

            GetCompositeScores(seqObj, charge, ms2ScanNum, out score, out nMatchedFragments);
            return new IcScores(nMatchedFragments, score, modifications);
        }

        public void GetCompositeScores(Sequence sequence, int parentIoncharge, int ms2ScanNum, out double score, out int nMatchedFragments)
        {
            score = 0d;
            nMatchedFragments = 0;

            var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return;

            var preFixIonCheck = new bool[sequence.Count + 1];
            var sufFixIonCheck = new bool[sequence.Count + 1];

            var scorer = new CompositeScorer(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, parentIoncharge+2), activationMethod: ActivationMethod);
            var cleavages = sequence.GetInternalCleavages();
            var cleavageIndex = 0;

            foreach (var c in cleavages)
            {
                bool prefixHit;
                bool suffixHit;
                score += scorer.GetFragmentScore(c.PrefixComposition, c.SuffixComposition, out prefixHit, out suffixHit);

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
    }

    public class IcScores
    {
        public IcScores(int nMatchedFragments, double score, string modifications)
        {
            NumMatchedFrags = nMatchedFragments;
            Score = score;
            Modifications = modifications;
        }

        public int NumMatchedFrags { get; private set; }
        public double Score { get; private set; } // this score is used to calculate p-value by generating function

        public string Modifications { get; private set; }

        public override string ToString()
        {
            return string.Join("\t",
                new[]
                {
                    NumMatchedFrags, Score,
                });
        }

        public static string GetScoreNames()
        {
            return "#MatchedFragments\tScore";
        }
    }
}
