using System;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class StatisticsTools
    {
        static public double GetLCCorrelation(Feature l, Feature r)
        {
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var llc = GetTruncatedLc(l, intersection);
            var rlc = GetTruncatedLc(r, intersection);
            return GetCorrelation(llc, rlc);
        }

        static public double GetIMSCorrelation(Feature l, Feature r)
        {
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var lims = GetTruncatedIms(l, intersection);
            var rims = GetTruncatedIms(r, intersection);
            return GetCorrelation(lims, rims);
        }


        static private double[] GetTruncatedLc(Feature feature, Rectangle intersection)
        {
            var sp = intersection.Left - feature.ScanLcStart;
            var ep = sp + intersection.Width;
            var lc = new double[intersection.Width];
            for (var i = sp; i < ep; i++)
                lc[i] = feature.LcApexPeakProfile[i];
            return lc;
        }

        static private double[] GetTruncatedIms(Feature feature, Rectangle intersection)
        {
            var sp = intersection.Top - feature.ScanImsStart;
            var ep = sp + intersection.Height;
            var ims = new double[intersection.Height];
            for (var i = sp; i < ep; i++)
                ims[i] = feature.ImsApexPeakProfile[i];
            return ims;
        }

        static public double GetIsotopeCorrelation(Feature[] f, double[] i)
        {
            var j = new double[f.Length];
            for (var k = 0; k < j.Length;k++ )
            {
                j[k] = f[k].IntensityMax;
            }
            return GetCorrelation(j, i);
        }
        
        static public double GetCorrelation(double[] v1, double[] v2)
        {
            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = Math.Sqrt(GetSampleVariance(v1, m1));
            var s2 = Math.Sqrt(GetSampleVariance(v2, m2));
            var rho = v1.Select((t, i) => (float)((t - m1) * (v2[i] - m2) / s1 / s2)).Sum();
            return Math.Max(0, rho / (v1.Length - 1));
        }

        static public double GetSampleMean(double[] x)
        {
            var m = x.Sum();
            return m / x.Length;
        }

        static public double GetSampleVariance(double[] x, double m)
        {
            var var = x.Sum(v => (v - m) * (v - m));
            return var / (x.Length - 1);
        }

    }
}
