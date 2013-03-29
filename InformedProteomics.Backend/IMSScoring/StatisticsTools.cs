using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class StatisticsTools
    {
        static public float GetLCCorrelation(Feature l, Feature r)
        {
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var llc = GetTruncatedLc(l, intersection);
            var rlc = GetTruncatedLc(r, intersection);
            return GetCorrelation(llc, rlc);
        }

        static public float GetIMSCorrelation(Feature l, Feature r)
        {
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var lims = GetTruncatedIms(l, intersection);
            var rims = GetTruncatedIms(r, intersection);
            return GetCorrelation(lims, rims);
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

        static public float GetCorrelation(float[] v1, float[] v2)
        {
            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = Math.Sqrt(GetSampleVariance(v1, m1));
            var s2 = Math.Sqrt(GetSampleVariance(v2, m2));
            var rho = v1.Select((t, i) => (float)((t - m1) * (v2[i] - m2) / s1 / s2)).Sum();
            return rho / (v1.Length - 1);
        }

        static public float GetSampleMean(float[] x)
        {
            var m = x.Sum();
            return m / x.Length;
        }

        static public float GetSampleVariance(float[] x, float m)
        {
            var var = x.Sum(v => (v - m) * (v - m));
            return var / (x.Length - 1);
        }
    }
}
