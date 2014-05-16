using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.BottomUp;
using InformedProteomics.Scoring.LikelihoodScoring;

namespace InformedProteomics.BottomUp.Scoring
{
    public class InformedBottomUpScorer
    {
        public InformedBottomUpScorer(LcMsRun run, AminoAcidSet aaSet, int minProductCharge, int maxProductCharge, Tolerance tolerance)
        {
            Run = run;
            AminoAcidSet = aaSet;
            MinProductCharge = minProductCharge;
            MaxProductCharge = maxProductCharge;
            Tolerance = tolerance;

            _rankScorer = new RankScore[4];
            for (var charge = 1; charge <= 4; charge++)
            {
                var scoringParamPath =
                    @"C:\cygwin\home\kims336\Data\QCShewQE\Scoring\HCD_QE_EDRN_RankProbabilities_Charge" + charge +
                    ".txt";
                _rankScorer[charge - 1] = new RankScore(scoringParamPath);
            }

            _scoredSpectra = new Dictionary<int, ScoredSpectrum>();
        }

        public LcMsRun Run { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinProductCharge { get; private set; }
        public int MaxProductCharge { get; private set; }
        public Tolerance Tolerance { get; private set; }

        private readonly RankScore[] _rankScorer;
        private readonly Dictionary<int, ScoredSpectrum> _scoredSpectra;

        public IcBottomUpScores GetScores(AminoAcid nTerm, string seqStr, AminoAcid cTerm, Composition composition, int charge, int ms2ScanNum)
        {
            ScoredSpectrum scoredSpectrum;
            var index = GetChargetScanNumPairIndex(charge, ms2ScanNum);
            if (!_scoredSpectra.TryGetValue(index, out scoredSpectrum))
            {
                var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (spec == null) return null;

                var rankScorer = charge > 4 ? _rankScorer[3] : _rankScorer[charge-1];
                scoredSpectrum = new ScoredSpectrum(spec, rankScorer, Tolerance);
                _scoredSpectra.Add(index, scoredSpectrum);
            }

            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, seqStr, AminoAcid.ProteinCTerm);
            if (seqGraph == null)
            {
                return null;
            }

            Tuple<double, string> scoreAndModifications = null;
            var bestScore = double.NegativeInfinity;
            var protCompositions = seqGraph.GetSequenceCompositions();
            for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);
                var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                if (!protCompositionWithH2O.Equals(composition)) continue;

                var curScoreAndModifications = seqGraph.GetScoreAndModifications(charge, scoredSpectrum);
                var curScore = curScoreAndModifications.Item1;
                if (curScore > bestScore)
                {
                    scoreAndModifications = curScoreAndModifications;
                    bestScore = curScore;
                }
            }

            if (scoreAndModifications == null) return null;

            var ms2Score = scoreAndModifications.Item1;
            var modifications = scoreAndModifications.Item2;

            return new IcBottomUpScores(ms2Score, modifications);
        }

        private int GetChargetScanNumPairIndex(int charge, int scanNum)
        {
            return charge * (Run.MaxLcScan + 1) + scanNum;
        }
    }
}
