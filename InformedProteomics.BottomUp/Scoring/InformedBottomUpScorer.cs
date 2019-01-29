using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.BottomUp;
using InformedProteomics.Scoring.LikelihoodScoring.Scoring;
using InformedProteomics.TopDown.Execution;

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
            // Ms2DetectorType.Orbitrap
            // Protocol.Standard
            _rankScorer = new RankScore(ActivationMethod.HCD, Enzyme.Trypsin);
            _scoredSpectra = new Dictionary<int, ScoredSpectrum>();
        }

        public LcMsRun Run { get; }
        public AminoAcidSet AminoAcidSet { get; }
        public int MinProductCharge { get; }
        public int MaxProductCharge { get; }
        public Tolerance Tolerance { get; }

        private readonly RankScore _rankScorer;
        private readonly Dictionary<int, ScoredSpectrum> _scoredSpectra;

        public IcBottomUpScores GetScores(DatabaseSequenceSpectrumMatch match, Composition composition, int charge,
            int ms2ScanNum)
        {
            return GetScores(match.Pre, match.Sequence, match.Post, match.NTerm, match.CTerm, composition, charge, ms2ScanNum);
        }

        public IcBottomUpScores GetScores(char pre, string sequence, char post, Composition composition, int charge,
            int ms2ScanNum)
        {
            var nTerm = pre == FastaDatabaseConstants.Delimiter ? AminoAcid.ProteinNTerm : AminoAcid.PeptideNTerm;
            var cTerm = post == FastaDatabaseConstants.Delimiter ? AminoAcid.ProteinCTerm : AminoAcid.PeptideCTerm;
            return GetScores(pre, sequence, post, nTerm, cTerm, composition, charge, ms2ScanNum);
        }

        public IcBottomUpScores GetScores(char pre, string sequence, char post, AminoAcid nTerm, AminoAcid cTerm, Composition composition, int charge, int ms2ScanNum)
        {
            var index = GetChargeScanNumPairIndex(charge, ms2ScanNum);
            if (!_scoredSpectra.TryGetValue(index, out var scoredSpectrum))
            {
                var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (spec == null) return null;
                scoredSpectrum = new ScoredSpectrum(spec, _rankScorer, charge, composition.Mass, Tolerance);
                _scoredSpectra.Add(index, scoredSpectrum);
            }

            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, nTerm, sequence, cTerm);
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

                var curScoreAndModifications = seqGraph.GetFragmentScoreAndModifications(scoredSpectrum);
                var curScore = curScoreAndModifications.Item1;
                if (curScore > bestScore)
                {
                    scoreAndModifications = curScoreAndModifications;
                    bestScore = curScore;
                }
            }

            if (scoreAndModifications == null) return null;

            var ms2Score = scoreAndModifications.Item1;

            // TODO: This assumes enzyme is trypsin
            const double probN = 0.99999;
            const double probC = 0.99999;
            const double sumAAProbabilities = 0.1;
            var creditN = Math.Log(probN / sumAAProbabilities);
            var penaltyN = Math.Log((1.0 - probN) / (1.0 - sumAAProbabilities));
            var creditC = Math.Log(probC / sumAAProbabilities);
            var penaltyC = Math.Log((1.0 - probC) / (1.0 - sumAAProbabilities));

            if (pre == 'K' || pre == 'R' || pre == FastaDatabaseConstants.Delimiter || pre == '-') ms2Score += creditN;
            else ms2Score += penaltyN;

            var lastResidue = sequence[sequence.Length - 1];
            if (lastResidue == 'K' || lastResidue == 'R' || post == FastaDatabaseConstants.Delimiter || post == '-') ms2Score += creditC;
            else ms2Score += penaltyC;

            var modifications = scoreAndModifications.Item2;

            return new IcBottomUpScores(ms2Score, modifications);
        }

        private int GetChargeScanNumPairIndex(int charge, int scanNum)
        {
            return charge * (Run.MaxLcScan + 1) + scanNum;
        }
    }
}
