using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.BottomUp;
using InformedProteomics.Scoring.LikelihoodScoring;
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
            const string scoringParamPath = @"\\protoapps\UserData\Sangtae\ScoringParams\HCD_QExactive_Tryp_RankProbabilities.txt";
            _rankScorer = new RankScore(scoringParamPath);
            _scoredSpectra = new Dictionary<int, ScoredSpectrum>();
        }

        public LcMsRun Run { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinProductCharge { get; private set; }
        public int MaxProductCharge { get; private set; }
        public Tolerance Tolerance { get; private set; }

        private readonly RankScore _rankScorer;
        private readonly Dictionary<int, ScoredSpectrum> _scoredSpectra;

        public IcBottomUpScores GetScores(SequenceSpectrumMatch match, Composition composition, int charge,
            int ms2ScanNum)
        {
            return GetScores(match.Pre, match.Sequence, match.Post, match.NTerm, match.CTerm, composition, charge, ms2ScanNum);
        }

        public IcBottomUpScores GetScores(char pre, string sequence, char post, Composition composition, int charge,
            int ms2ScanNum)
        {
            var nTerm = pre == FastaDatabase.Delimiter ? AminoAcid.ProteinNTerm : AminoAcid.PeptideNTerm;
            var cTerm = post == FastaDatabase.Delimiter ? AminoAcid.ProteinCTerm : AminoAcid.PeptideCTerm;
            return GetScores(pre, sequence, post, nTerm, cTerm, composition, charge, ms2ScanNum);
        }

        public IcBottomUpScores GetScores(char pre, string sequence, char post, AminoAcid nTerm, AminoAcid cTerm, Composition composition, int charge, int ms2ScanNum)
        {
            ScoredSpectrum scoredSpectrum;
            var index = GetChargetScanNumPairIndex(charge, ms2ScanNum);
            if (!_scoredSpectra.TryGetValue(index, out scoredSpectrum))
            {
                var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (spec == null) return null;
                scoredSpectrum = new ScoredSpectrum(spec, _rankScorer, charge, Tolerance);
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

            // TODO: This assumes enzyme is trypsin
            if (pre == 'K' || pre == 'R' || pre == FastaDatabase.Delimiter || pre == '-') ms2Score += 2.21;
            else ms2Score += -0.14;

            var lastResidue = sequence[sequence.Length - 1];
            if (lastResidue == 'K' || lastResidue == 'R' || post == FastaDatabase.Delimiter || post == '-') ms2Score += 2.29;
            else ms2Score += -2.41;

            var modifications = scoreAndModifications.Item2;

            return new IcBottomUpScores(ms2Score, modifications);
        }

        private int GetChargetScanNumPairIndex(int charge, int scanNum)
        {
            return charge * (Run.MaxLcScan + 1) + scanNum;
        }
    }
}
