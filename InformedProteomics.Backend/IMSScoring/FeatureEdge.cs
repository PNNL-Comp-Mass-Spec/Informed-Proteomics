using System;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureEdge
    {
        public FeatureNode LNode { get; private set; }
        public FeatureNode RNode { get; private set; }
        public float Weight { get; private set; } // used to calculate weight of a path
        public float Score { get; private set; } // used to calculate score
        private readonly int _ratio;
        private readonly float _lcCorrelation;
        private readonly float _imsCorrelation;

        public FeatureEdge(FeatureNode l, FeatureNode r)
        {
            LNode = l;
            RNode = r;

            _ratio = GetRatio(l.Feature.IntensityMax, r.Feature.IntensityMax);
            _lcCorrelation = GetLCCorrelation(l, r);
            _imsCorrelation = GetIMSCorrelation(l, r);
           
            Weight = GetWeight();
            Score = GetScore();
        }

        static private float GetLCCorrelation(FeatureNode l, FeatureNode r)
        {
            var intersection = Rectangle.Intersect(l.Feature.GetBoundary(), r.Feature.GetBoundary());
            var llc = GetTruncatedLc(l.Feature, intersection);
            var rlc = GetTruncatedLc(r.Feature, intersection);
            return GetCorrelation(llc, rlc);
        }

        static private float GetIMSCorrelation(FeatureNode l, FeatureNode r)
        {
            var intersection = Rectangle.Intersect(l.Feature.GetBoundary(), r.Feature.GetBoundary());
            var lims = GetTruncatedIms(l.Feature, intersection);
            var rims = GetTruncatedIms(r.Feature, intersection);
            return GetCorrelation(lims, rims);
        }

        private float GetWeight()
        {
            return SubScoreFactory.GetKLDivergence(LNode.FragmentIonClassBase, RNode.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, LNode.Parameter);
        }

        private float GetScore()
        {
            var score = SubScoreFactory.GetRatioScore(LNode.FragmentIonClassBase, RNode.FragmentIonClassBase, _ratio, LNode.Parameter);
            score += SubScoreFactory.GetLCCorrelationScore(LNode.FragmentIonClassBase, RNode.FragmentIonClassBase, _lcCorrelation, LNode.Parameter);
            score += SubScoreFactory.GetIMSCorrelationScore(LNode.FragmentIonClassBase, RNode.FragmentIonClassBase, _imsCorrelation, LNode.Parameter);

            return score;
        }

        static private float[] GetTruncatedLc(Feature feature, Rectangle intersection)
        {
            var sp = intersection.Left - feature.ScanLcStart;
            var ep = sp + intersection.Width;
            var lc = new float[intersection.Width];
            for (var i = sp; i < ep; i++)
                lc[i] = feature.LcApexPeakProfile[i];
            return lc;
        }

        static private float[] GetTruncatedIms(Feature feature, Rectangle intersection)
        {
            var sp = intersection.Top - feature.ScanImsStart;
            var ep = sp + intersection.Height;
            var ims = new float[intersection.Height];
            for (var i = sp; i < ep; i++)
                ims[i] = feature.ImsApexPeakProfile[i];
            return ims;
        }

        static private int GetRatio(double v1, double v2)
        {
            if (v1 <= 0)
            {
                if (v2 <= 0) return -11;
                return -12;
            }
            if (v2 <= 0) return 11;
            double r;
            var f = 1;
            if (v1 > v2) r = v1 / v2;
            else
            {
                r = v2 / v1;
                f = -1;
            }
            return (int)(Math.Min(10, r) * f);
        }

        static private float GetCorrelation(float[] v1, float[] v2)
        {
            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = Math.Sqrt(GetSampleVariance(v1, m1));
            var s2 = Math.Sqrt(GetSampleVariance(v2, m2));
            var rho = v1.Select((t, i) => (float) ((t - m1)*(v2[i] - m2)/s1/s2)).Sum();
            return rho/(v1.Length - 1);
        }

        static private float GetSampleMean(float[] x)
        {
            var m = x.Sum();
            return m / x.Length;
        }

        static private float GetSampleVariance(float[] x, float m)
        {
            var var = x.Sum(v => (v - m) * (v - m));
            return var / (x.Length - 1);
        }
    }
}
