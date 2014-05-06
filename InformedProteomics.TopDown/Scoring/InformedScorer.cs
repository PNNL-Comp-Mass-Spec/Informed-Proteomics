using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class InformedScorer
    {
        public InformedScorer(LcMsRun run, AminoAcidSet aaSet, int minProductCharge, int maxProductCharge, Tolerance tolerance, double ms2CorrThreshold = 0.7)
        {
            Run = run;
            AminoAcidSet = aaSet;
            MinProductCharge = minProductCharge;
            MaxProductCharge = maxProductCharge;
            Tolerance = tolerance;
            Ms2CorrThreshold = ms2CorrThreshold;
        }

        public LcMsRun Run { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinProductCharge { get; private set; }
        public int MaxProductCharge { get; private set; }
        public Tolerance Tolerance { get; private set; }
        public double Ms2CorrThreshold { get; private set; }

        public IcScores GetScores(string seqStr, Composition composition, int charge, int ms2ScanNum)
        {
            var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return null;

            var annotation = "_." + seqStr + "._";
            var scorer = new CorrMatchedPeakCounter(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, charge), Ms2CorrThreshold);

            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, annotation);
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
                
                var curScoreAndModifications = seqGraph.GetScoreAndModifications(charge, scorer);
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

            var ion = new Ion(composition, charge);

            // IsotopeCorrPrevMs1
            var precursorSpec = Run.GetSpectrum(Run.GetPrecursorScanNum(ms2ScanNum));
            var isotopeCorrPrevMs1 = precursorSpec == null ? 0.0 : precursorSpec.GetCorrScore(ion, Tolerance);

            // IsotopeCorrNextMs1
            var nextMs1Spec = Run.GetSpectrum(Run.GetNextScanNum(ms2ScanNum, 1));
            var isotopeCorrNextMs1 = nextMs1Spec == null ? 0.0 : nextMs1Spec.GetCorrScore(ion, Tolerance);

            var mostAbundantIsotopeMz = ion.GetMostAbundantIsotopeMz();
            var xicThisPeak = Run.GetExtractedIonChromatogram(mostAbundantIsotopeMz, Tolerance, ms2ScanNum);
            if (xicThisPeak.Count < 2)
            {
                return new IcScores(ms2Score, isotopeCorrPrevMs1, isotopeCorrNextMs1, 0.0, 0.0, 0.0, modifications);
            }

            // check whether next isotope peak exists
            var nextIsotopeMz = mostAbundantIsotopeMz + Constants.C13MinusC12 / charge;
            var xicNextIsotope = Run.GetExtractedIonChromatogram(nextIsotopeMz, Tolerance, ms2ScanNum);
            var corrMostAbundantPlusOneIsotope = xicThisPeak.GetCorrelation(xicNextIsotope);

            var mzChargeMinusOne = new Ion(composition, charge - 1).GetMostAbundantIsotopeMz();
            var xicMinusOneCharge = Run.GetExtractedIonChromatogram(mzChargeMinusOne, Tolerance, ms2ScanNum);
            var corrMinusOneCharge = xicMinusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicMinusOneCharge) : 0;

            var mzChargePlusOne = new Ion(composition, charge + 1).GetMostAbundantIsotopeMz();
            var xicPlusOneCharge = Run.GetExtractedIonChromatogram(mzChargePlusOne, Tolerance, ms2ScanNum);
            var corrPlusOneCharge = xicPlusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicPlusOneCharge) : 0;

            return new IcScores(ms2Score, isotopeCorrPrevMs1, isotopeCorrNextMs1, corrMostAbundantPlusOneIsotope, corrMinusOneCharge, corrPlusOneCharge, modifications);
        }
    }

    public class IcScores
    {
        public IcScores(
            double ms2Score, 
            double isotopeCorrPrevMs1, double isotopeCorrNextMs1, 
            double corrMostAbundantPlusOneIsotope, 
            double chargeCorrMinusOne, double chargeCorrPlusOne,
            string modifications)
        {
            Ms2Score = ms2Score;
            IsotopeCorrPrevMs1 = isotopeCorrPrevMs1;
            IsotopeCorrNextMs1 = isotopeCorrNextMs1;
            CorrMostAbundantPlusOneIsotope = corrMostAbundantPlusOneIsotope;
            ChargeCorrMinusOne = chargeCorrMinusOne;
            ChargeCorrPlusOne = chargeCorrPlusOne;
            Modifications = modifications;
        }

        public double Ms2Score { get; private set; }
        public double IsotopeCorrPrevMs1 { get; private set; }
        public double IsotopeCorrNextMs1 { get; private set; }
        public double CorrMostAbundantPlusOneIsotope { get; private set; }
        public double ChargeCorrMinusOne { get; private set; }
        public double ChargeCorrPlusOne { get; private set; }
        public string Modifications { get; private set; }

        public override string ToString()
        {
            return string.Join("\t",
                new[]
                {
                    Ms2Score, IsotopeCorrPrevMs1, IsotopeCorrNextMs1, CorrMostAbundantPlusOneIsotope, ChargeCorrMinusOne,
                    ChargeCorrPlusOne
                });
        }

        public static string GetScoreNames()
        {
            return
                "Ms2Score\tIsotopeCorrPrevMs1\tIsotopeCorrNextMs1\tCorrMostAbundantPlusOneIsotope\tChargeCorrMinusOne\tChargeCorrPlusOne";
        }
    }

    //public IcScores GetScores(Sequence sequence, int charge, int ms2ScanNum)
    //{
    //    var spec = Run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
    //    if (spec == null) return null;

    //    var scorer = new CorrMatchedPeakCounter(spec, Tolerance, MinProductCharge, Math.Min(MaxProductCharge, charge), Ms2CorrThreshold);
    //    var internalCleavages = sequence.GetInternalCleavages();
    //    var ms2Score = internalCleavages.Sum(
    //        cleavage => scorer.GetFragmentScore(cleavage.PrefixComposition, cleavage.SuffixComposition));

    //    var sequenceCompositionWithH2O = sequence.Composition + Composition.H2O;
    //    var ion = new Ion(sequenceCompositionWithH2O, charge);

    //    // IsotopeCorrPrevMs1
    //    var precursorSpec = Run.GetSpectrum(Run.GetPrecursorScanNum(ms2ScanNum));
    //    var isotopeCorrPrevMs1 = precursorSpec == null ? 0.0 : precursorSpec.GetCorrScore(ion, Tolerance);

    //    // IsotopeCorrNextMs1
    //    var nextMs1Spec = Run.GetSpectrum(Run.GetNextScanNum(ms2ScanNum, 1));
    //    var isotopeCorrNextMs1 = nextMs1Spec == null ? 0.0 : nextMs1Spec.GetCorrScore(ion, Tolerance);

    //    var mostAbundantIsotopeMz = ion.GetMostAbundantIsotopeMz();
    //    var xicThisPeak = Run.GetExtractedIonChromatogram(mostAbundantIsotopeMz, Tolerance, ms2ScanNum);
    //    if (xicThisPeak.Count < 2)
    //    {
    //        return new IcScores(ms2Score, isotopeCorrPrevMs1, isotopeCorrNextMs1, 0.0, 0.0, 0.0);
    //    }

    //    // check whether next isotope peak exists
    //    var nextIsotopeMz = mostAbundantIsotopeMz + Constants.C13MinusC12 / charge;
    //    var xicNextIsotope = Run.GetExtractedIonChromatogram(nextIsotopeMz, Tolerance, ms2ScanNum);
    //    var corrMostAbundantPlusOneIsotope = xicThisPeak.GetCorrelation(xicNextIsotope);

    //    var mzChargeMinusOne = new Ion(sequenceCompositionWithH2O, charge - 1).GetMostAbundantIsotopeMz();
    //    var xicMinusOneCharge = Run.GetExtractedIonChromatogram(mzChargeMinusOne, Tolerance, ms2ScanNum);
    //    var corrMinusOneCharge = xicMinusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicMinusOneCharge) : 0;

    //    var mzChargePlusOne = new Ion(sequenceCompositionWithH2O, charge + 1).GetMostAbundantIsotopeMz();
    //    var xicPlusOneCharge = Run.GetExtractedIonChromatogram(mzChargePlusOne, Tolerance, ms2ScanNum);
    //    var corrPlusOneCharge = xicPlusOneCharge.Count >= 3 ? xicThisPeak.GetCorrelation(xicPlusOneCharge) : 0;

    //    return new IcScores(ms2Score, isotopeCorrPrevMs1, isotopeCorrNextMs1, corrMostAbundantPlusOneIsotope, corrMinusOneCharge, corrPlusOneCharge);
    //}
}
