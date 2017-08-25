using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MathAndStats;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using InformedProteomics.Backend.Data.Composition;

    public abstract class AbstractFragmentScorer : IScorer
    {
        protected AbstractFragmentScorer(Spectrum spec, Tolerance tol, int minCharge = 1, int maxCharge = 20, double relativeIsotopeIntensityThreshold = 0.7, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            Ms2Spectrum = spec;
            Tolerance = tol;
            MinProductCharge = minCharge;
            MaxProductCharge = maxCharge;

            ActivationMethod actToUse = activationMethod;
            if (actToUse == ActivationMethod.Unknown)
            {
                var spectrum = spec as ProductSpectrum;
                if (spectrum != null) actToUse = spectrum.ActivationMethod;
                else actToUse = ActivationMethod.CID;
            }

            if (actToUse == ActivationMethod.CID || actToUse == ActivationMethod.HCD) BaseIonTypes = BaseIonTypesCID;
            else if (actToUse == ActivationMethod.ETD) BaseIonTypes = BaseIonTypesETD;
            else if (actToUse == ActivationMethod.UVPD) BaseIonTypes = BaseIonTypesUVPD;
            else BaseIonTypes = BaseIonTypesCID;

            if (actToUse == ActivationMethod.UVPD)
            {
                this.PrefixOffsetMass = this.BaseIonTypes[0].OffsetComposition.Mass;
                this.SuffixOffsetMass = this.BaseIonTypes[2].OffsetComposition.Mass;
            }
            else
            {
                this.PrefixOffsetMass = this.BaseIonTypes[0].OffsetComposition.Mass;
                this.SuffixOffsetMass = this.BaseIonTypes[1].OffsetComposition.Mass;
            }

            RelativeIsotopeIntensityThreshold = relativeIsotopeIntensityThreshold;
        }

        public abstract double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null);

        protected IntRange GetMinMaxChargeRange(Composition fragmentComposition)
        {
            var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();
            var fragmentIonMostAbuMass = fragmentComposition.Mass + Constants.C13MinusC12 * mostAbundantIsotopeIndex;
            return GetMinMaxChargeRange(fragmentIonMostAbuMass);
        }

        protected IntRange GetMinMaxChargeRange(double fragmentIonMostAbuMass)
        {
            var minMz = Ms2Spectrum.Peaks.First().Mz;
            var maxMz = Ms2Spectrum.Peaks.Last().Mz;

            var maxCharge = (int)Math.Floor(fragmentIonMostAbuMass / (minMz - Constants.Proton));
            var minCharge = (int)Math.Ceiling(fragmentIonMostAbuMass / (maxMz - Constants.Proton));
            if (maxCharge < 1 || maxCharge > MaxProductCharge) maxCharge = MaxProductCharge;
            if (minCharge < 1 || minCharge < MinProductCharge) minCharge = MinProductCharge;

            return new IntRange(minCharge, maxCharge);
        }

        protected IEnumerable<DeconvolutedPeak> FindMatchedPeaks(Composition fragmentComposition,
                   double corrThreshold, double distThreshold)
        {
            var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();
            var fragmentIonMass = fragmentComposition.Mass;

            //var matchedPeak = new MatchedFragmentPeak();
            //var deconvPeak = new DeconvolutedPeak()
            if (fragmentIonMass < Ms2Spectrum.Peaks.First().Mz) yield break;

            var prevObservedCharge = 0;
            var fragmentIonMostAbuMass = fragmentIonMass + Constants.C13MinusC12 * mostAbundantIsotopeIndex;
            var chargeRange = GetMinMaxChargeRange(fragmentIonMostAbuMass);

            for (var charge = chargeRange.Min; charge <= chargeRange.Max; charge++)
            {
                var ion = new Ion(fragmentComposition, charge);

                var observedPeaks = Ms2Spectrum.GetAllIsotopePeaks(ion, Tolerance, RelativeIsotopeIntensityThreshold);
                if (observedPeaks == null)
                {
                    if (prevObservedCharge > 0 && charge - prevObservedCharge > 1) yield break;
                    continue;
                }

                var distCorr = GetDistCorr(ion, observedPeaks);
                if (distCorr.Item2 < corrThreshold && distCorr.Item1 > distThreshold)
                {
                    if (prevObservedCharge > 0 && charge - prevObservedCharge > 1) yield break;
                    continue;
                }
                var matchedPeak = new DeconvolutedPeak(fragmentIonMass, observedPeaks[mostAbundantIsotopeIndex].Intensity, charge, distCorr.Item2, distCorr.Item1, observedPeaks);
                prevObservedCharge = charge;

                yield return matchedPeak;
            }
        }
        /*
        protected Peak[] FindMostAbundantPeak(Composition.Composition fragmentComposition,
                   double corrThreshold, double distThreshold,
                   out int observedCharge, out double envelopeCorr, out double envelopeDist)
        {
            //Peak[] intenseObservedPeaks = null;
            var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();
            var fragmentIonMass = fragmentComposition.Mass;
            observedCharge = 0;
            envelopeCorr = 0d;
            envelopeDist = 1.0d;

            if (fragmentIonMass < Ms2Spectrum.Peaks.First().Mz) return null;

            var fragmentIonMostAbuMass = fragmentIonMass + Constants.C13MinusC12 * mostAbundantIsotopeIndex;
            var chargeRange = GetMinMaxChargeRange(fragmentIonMostAbuMass);

            for (var charge = chargeRange.Min; charge <= chargeRange.Max; charge++)
            {
                var ion = new Ion(fragmentComposition, charge);

                var observedPeaks = Ms2Spectrum.GetAllIsotopePeaks(ion, Tolerance, RelativeIsotopeIntensityThreshold);
                if (observedPeaks == null) continue;

                var distCorr = GetDistCorr(ion, observedPeaks);
                if (distCorr.Item2 < corrThreshold && distCorr.Item1 > distThreshold) continue;
                //var mostAbuPeak = observedPeaks[mostAbundantIsotopeIndex];
                //if (intenseObservedPeaks == null || mostAbuPeak.Intensity > intenseObservedPeaks[mostAbundantIsotopeIndex].Intensity)
                //{
                    //intenseObservedPeaks = observedPeaks;
                observedCharge = charge;
                envelopeDist = distCorr.Item1;
                envelopeCorr = distCorr.Item2;
                return observedPeaks;
                //}
            }
            return null;
        }
        */
        protected Peak[] FindMostIntensePeak(Composition fragmentComposition, double corrThreshold, double distThreshold, out int observedCharge, out double envelopeCorr, out double envelopeDist)
        {
            Peak[] intenseObservedPeaks = null;
            var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();
            var fragmentIonMass = fragmentComposition.Mass;
            observedCharge = 0;
            envelopeCorr = 0d;
            envelopeDist = 1.0d;

            if (fragmentIonMass < Ms2Spectrum.Peaks.First().Mz) return null;

            var fragmentIonMostAbuMass = fragmentIonMass + Constants.C13MinusC12 * mostAbundantIsotopeIndex;
            var chargeRange = GetMinMaxChargeRange(fragmentIonMostAbuMass);

            for (var charge = chargeRange.Min; charge <= chargeRange.Max; charge++)
            {
                var ion = new Ion(fragmentComposition, charge);

                var observedPeaks = Ms2Spectrum.GetAllIsotopePeaks(ion, Tolerance, RelativeIsotopeIntensityThreshold);
                if (observedPeaks == null) continue;

                var distCorr = GetDistCorr(ion, observedPeaks);
                if (distCorr.Item2 < corrThreshold && distCorr.Item1 > distThreshold) continue;
                var mostAbuPeak = observedPeaks[mostAbundantIsotopeIndex];

                if (intenseObservedPeaks == null || mostAbuPeak.Intensity > intenseObservedPeaks[mostAbundantIsotopeIndex].Intensity)
                {
                    intenseObservedPeaks = observedPeaks;
                    observedCharge = charge;
                    envelopeDist = distCorr.Item1;
                    envelopeCorr = distCorr.Item2;
                }
            }
            return intenseObservedPeaks;
        }

        public static Tuple<double, double> GetDistCorr(Ion ion, Peak[] observedPeaks)
        {
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
            var envelope = isotopomerEnvelope;
            var observedIntensities = new double[envelope.Length];
            for (var i = 0; i < isotopomerEnvelope.Length; i++)
            {
                if (observedPeaks[i] != null) observedIntensities[i] = observedPeaks[i].Intensity;
            }

            return FitScoreCalculator.GetDistanceAndCorrelation(envelope, observedIntensities);
            //var bcDist = FitScoreCalculator.GetBhattacharyyaDistance(envelope, observedIntensities);
            //var corr = FitScoreCalculator.GetPearsonCorrelation(envelope, observedIntensities);
            //return new Tuple<double, double>(bcDist, corr);
        }

        public readonly Spectrum Ms2Spectrum;
        protected readonly Tolerance Tolerance;
        protected readonly int MinProductCharge;
        protected readonly int MaxProductCharge;
        protected readonly BaseIonType[] BaseIonTypes;

        protected double PrefixOffsetMass;

        protected double SuffixOffsetMass;

        protected double RelativeIsotopeIntensityThreshold;

        protected static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD, BaseIonTypesUVPD;
        static AbstractFragmentScorer()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
            BaseIonTypesUVPD = new[] { BaseIonType.A, BaseIonType.Ar, BaseIonType.B, BaseIonType.C, BaseIonType.X, BaseIonType.Xr, BaseIonType.Y, BaseIonType.YM1, BaseIonType.Z, BaseIonType.Zr };
        }
    }
}
