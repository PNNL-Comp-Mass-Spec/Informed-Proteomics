using System;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Financial;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.TopDown.Scoring
{
    using InformedProteomics.Backend.Data.Sequence;

    public class LikelihoodScorer : AbstractFragmentScorer
    {
        public const double CorrThreshold = 0.7;
        public const double DistThreshold = 0.07;
        private readonly bool _includeMassErrorScore;

        public LikelihoodScorer(LikelihoodScoringModel model, ProductSpectrum ms2Spec, Tolerance tolerance, int minCharge, int maxCharge, bool massErrorScore = true)
            : base(ms2Spec, tolerance, minCharge, maxCharge, 0.1)
        {
            //_model = model;

            //var refIntensity = ms2Spec.Peaks.Max(p => p.Intensity) * 0.1;
            //var medIntensity = ms2Spec.Peaks.Select(p => p.Intensity).Median();
            //_refIntensity = Math.Min(medIntensity, refIntensity);

            _refIntensity = GetRefIntensity(ms2Spec.Peaks);

            _includeMassErrorScore = massErrorScore;
        }

        public static double GetRefIntensity(Peak[] peaks)
        {
            return peaks.Sum(p => p.Intensity);
        }

        public override double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
        {
            var score = 0.0;
            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;

                var observedCharge = 0;
                var envelopeCorr = 0d;
                var envelopeDist = 0d;
                var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();

                var observedPeaks = FindMostIntensePeak(fragmentComposition, CorrThreshold, DistThreshold, out observedCharge,
                    out envelopeCorr, out envelopeDist);
                var fragmentIonMass = fragmentComposition.Mass;

                if (observedPeaks == null) continue;
                var observedMostAbuPeak = observedPeaks[mostAbundantIsotopeIndex];

                var observedMass = Ion.GetMonoIsotopicMass(observedMostAbuPeak.Mz, observedCharge, mostAbundantIsotopeIndex);
                var massErrorPpm = (Math.Abs(observedMass - fragmentIonMass)/fragmentIonMass)*1e6;

                score += 1;
                var intScore = (observedMostAbuPeak.Intensity / _refIntensity) * 10;
                var corrScore = (fragmentIonMass > 1300 & envelopeCorr > 0.7) ? (envelopeCorr - 0.7) : 0;
                var distScore = (fragmentIonMass > 1300 & envelopeDist < 0.07) ? 0.3 - 3.75 * envelopeDist : 0;
                score += intScore;
                score += corrScore;
                score += distScore;

                if (_includeMassErrorScore)
                {
                    /*score += _model.GetNodeScore(Ms2Spectrum.ActivationMethod, baseIonType.IsPrefix,
                        fragmentIonMass, observedCharge,
                        envelopeCorr, envelopeDist,
                        observedMostAbuPeak.Intensity / _refIntensity, massErrorPpm);
                     */
                    var massErrorScore = GetMassErrorScore(massErrorPpm);
                    score += massErrorScore;
                }
                else
                {
                    /*score += _model.GetNodeScoreWithoutMassError(Ms2Spectrum.ActivationMethod, baseIonType.IsPrefix,
                        fragmentIonMass, observedCharge,
                        envelopeCorr, envelopeDist,
                        observedMostAbuPeak.Intensity / _refIntensity);
                     */
                }
                //score += _model.GetScore(baseIonType, bestCorrScore, bestObsIntensity);
            }
            return score;
        }

        public double GetMassErrorScore(double massErrorPpm)
        {
            //y = normpdf(bin_Z, 5.0140, 3.1534);
            //y2 = normpdf(bin_Z, 6.3361, 4.0447);
            var p = Normal.PDF(5.014, 3.1534, massErrorPpm);
            var pnull = Normal.PDF(7.3361, 4.0447, massErrorPpm);
            return Math.Max(Math.Log(p / pnull), 0);
        }

        public Tuple<double, int> GetLikelihoodAndPeakCountScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
        {
            var score = 0.0;
            var peakCount = 0;
            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;

                var observedCharge = 0;
                var envelopeCorr = 0d;
                var envelopeDist = 0d;
                var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();

                var observedPeaks = FindMostIntensePeak(fragmentComposition, CorrThreshold, DistThreshold, out observedCharge,
                    out envelopeCorr, out envelopeDist);
                var fragmentIonMass = fragmentComposition.Mass;

                if (observedPeaks == null) continue;
                var observedMostAbuPeak = observedPeaks[mostAbundantIsotopeIndex];

                var observedMass = Ion.GetMonoIsotopicMass(observedMostAbuPeak.Mz, observedCharge, mostAbundantIsotopeIndex);
                var massErrorPpm = (Math.Abs(observedMass - fragmentIonMass) / fragmentIonMass) * 1e6;

                peakCount++;

                score += 1;
                var intScore = (observedMostAbuPeak.Intensity / _refIntensity) * 10;
                var corrScore = (fragmentIonMass > 1300 & envelopeCorr > 0.7) ? (envelopeCorr - 0.7) : 0;
                var distScore = (fragmentIonMass > 1300 & envelopeDist < 0.07) ? 0.3 - 3.75 * envelopeDist : 0;
                score += intScore;
                score += corrScore;
                score += distScore;
                if (_includeMassErrorScore)
                {
                    var massErrorScore = GetMassErrorScore(massErrorPpm);
                    score += massErrorScore;
                }

                /*
                if (_includeMassErrorScore)
                {
                    score += _model.GetNodeScore(Ms2Spectrum.ActivationMethod, baseIonType.IsPrefix,
                        fragmentIonMass, observedCharge,
                        envelopeCorr, envelopeDist,
                        observedMostAbuPeak.Intensity / _refIntensity, massErrorPpm);
                }
                else
                {
                    score += _model.GetNodeScoreWithoutMassError(Ms2Spectrum.ActivationMethod, baseIonType.IsPrefix,
                        fragmentIonMass, observedCharge,
                        envelopeCorr, envelopeDist,
                        observedMostAbuPeak.Intensity / _refIntensity);
                }*/

                //score += _model.GetScore(baseIonType, bestCorrScore, bestObsIntensity);
            }
            return new Tuple<double, int>(score, peakCount);
        }

        //private readonly LikelihoodScoringModel _model;
        private readonly double _refIntensity;
    }
}
