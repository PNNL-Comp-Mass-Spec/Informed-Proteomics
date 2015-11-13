using System;
using System.CodeDom;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring
{
    public class CompositeScorer : AbstractFragmentScorer
    {
        public CompositeScorer(Spectrum ms2Spec, Tolerance tol, int minCharge, int maxCharge, double relativeIsotopeIntensityThreshold = 0.1)
            : base(ms2Spec, tol, minCharge, maxCharge, relativeIsotopeIntensityThreshold)
        {
            ReferencePeakIntensity = GetRefIntensity(ms2Spec.Peaks);
        }

        public CompositeScorer(Spectrum ms2Spec, Tolerance tol, double relativeIsotopeIntensityThreshold = 0.1)
            : base(ms2Spec, tol, 1, 20, relativeIsotopeIntensityThreshold)
        {
        }

        public const double CorrThreshold = 0.7;
        public const double DistThreshold = 0.03;
        public double ReferencePeakIntensity { get; protected set; }

        public override double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition)
        {
            bool prefixHit, suffixHit;
            return GetFragmentScore(prefixFragmentComposition, suffixFragmentComposition, out prefixHit, out suffixHit);
        }
        
        public double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition, out bool prefixHit, out bool suffixHit)
        {
            var score = 0.0;
            
            prefixHit = false;
            suffixHit = false;

            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                              ? prefixFragmentComposition + baseIonType.OffsetComposition
                              : suffixFragmentComposition + baseIonType.OffsetComposition;

                if (fragmentComposition.Mass < Ms2Spectrum.Peaks[0].Mz) continue;

                var param = baseIonType.IsPrefix ? ScoreParam.Prefix : ScoreParam.Suffix;
                var fragmentIonMass = fragmentComposition.Mass;
                
                var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();

                foreach (var matchedPeak in FindMatchedPeaks(fragmentComposition, CorrThreshold, DistThreshold))
                {
                    var observedMostAbuPeak = matchedPeak.ObservedPeaks[mostAbundantIsotopeIndex];
                    var observedMass = Ion.GetMonoIsotopicMass(observedMostAbuPeak.Mz, matchedPeak.Charge, mostAbundantIsotopeIndex);
                    var massErrorPpm = (Math.Abs(observedMass - fragmentIonMass) / fragmentIonMass) * 1e6;
                    
                    score += param.Count;
                    score += param.Intensity * Math.Min(observedMostAbuPeak.Intensity / ReferencePeakIntensity, 1.0); // intensity-based scoring
                    score += param.Dist * matchedPeak.Dist; // Envelope distance-based scoring
                    score += param.Corr * matchedPeak.Corr; // Envelope correlation-based scoring
                    score += param.MassError * massErrorPpm; // Envelope correlation-based scoring

                    if (baseIonType.IsPrefix) prefixHit = true;
                    else suffixHit = true;
                }
            }

            if (prefixHit && suffixHit) score += ScoreParam.ComplementaryIonCount;
            return score;
        }

        public static double GetRefIntensity(Peak[] peaks)
        {
            var intensities = peaks.Select(p => p.Intensity).ToArray();
            Array.Sort(intensities, 0, intensities.Length);
            var refIntensity = intensities[(int)(0.8750 * intensities.Length)];
            return refIntensity;
        }

        public static double GetProbability(double score)
        {
            var eta = Math.Exp(ScoreParam.Beta0 + score);
            return eta/(eta + 1);
        }

        internal static ScoreWeight ScoreParam; // score weights without mass error temrs for generating function evaluation
        private const double WeightScaleFactor = 4.0;

        static CompositeScorer()
        {
            /*ScoreParam = new ScoreWeight()
            {
                Beta0 = -5.39335774417474*WeightScaleFactor,
                Cutoff = 2.5 * WeightScaleFactor,
                ComplementaryIonCount = -0.378768218613132 * WeightScaleFactor,
                Prefix = new IonScoreWeight
                {
                    Count = -1.125486077 * WeightScaleFactor,
                    ConsecutiveMatch = 0.125523007 * WeightScaleFactor,
                    Intensity = 0.323923971 * WeightScaleFactor,
                    Corr = 1.473096149 * WeightScaleFactor,
                    Dist = 0.179502658 * WeightScaleFactor,
                    MassError = -0.032959575 * WeightScaleFactor,
                },
                Suffix = new IonScoreWeight()
                {
                    Count = -0.886170429 * WeightScaleFactor,
                    ConsecutiveMatch = 0.083704707 * WeightScaleFactor,
                    Intensity = 0.372673677 * WeightScaleFactor,
                    Corr = 1.18929338 * WeightScaleFactor,
                    Dist = -0.284321389 * WeightScaleFactor,
                    MassError = -0.026141488 * WeightScaleFactor,
                },
            };*/
            var trainedParam = new double[]
            {
                -4.35615875783574,
                -0.263297798433790,
                -0.863126458013878,
                0.238325241230975,
                0.294277332664349,
                0.882288576167220,
                0.732263929496114,
                -0.0199416790747384,
                -0.498811136842684,
                0.180874679890120,
                0.299790673329884,
                0.492391631834416,
                0.314410358803801,
                -0.0134845040032297,
            };

            ScoreParam = new ScoreWeight()
            {
                Beta0 = trainedParam[0] * WeightScaleFactor,
                Cutoff = 2.0 * WeightScaleFactor,
                ComplementaryIonCount = trainedParam[1] * WeightScaleFactor,
                Prefix = new IonScoreWeight
                {
                    Count = trainedParam[2] * WeightScaleFactor,
                    ConsecutiveMatch = trainedParam[3] * WeightScaleFactor,
                    Intensity = trainedParam[4] * WeightScaleFactor,
                    Corr = trainedParam[5] * WeightScaleFactor,
                    Dist = trainedParam[6] * WeightScaleFactor,
                    MassError = trainedParam[7] * WeightScaleFactor,
                },
                Suffix = new IonScoreWeight()
                {
                    Count = trainedParam[8] * WeightScaleFactor,
                    ConsecutiveMatch = trainedParam[9] * WeightScaleFactor,
                    Intensity = trainedParam[10] * WeightScaleFactor,
                    Corr = trainedParam[11] * WeightScaleFactor,
                    Dist = trainedParam[12] * WeightScaleFactor,
                    MassError = trainedParam[13] * WeightScaleFactor,
                },
            };
        }

        internal class ScoreWeight
        {
            internal double ComplementaryIonCount;
            internal IonScoreWeight Prefix;
            internal IonScoreWeight Suffix;
            internal double Cutoff;

            internal double Beta0;
        }

        internal class IonScoreWeight
        {
            internal double Count;
            
            internal double Intensity;
            internal double Corr;
            internal double Dist;
            internal double MassError;

            internal double ConsecutiveMatch;
        }
    }
}
