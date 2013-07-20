using System;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class StatisticsTools
    {
        static public double GetLcCorrelation(Feature l, Feature r)
        {
            if (l == null || r == null) return -1;
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var llc = GetTruncatedLc(l, intersection);
            var rlc = GetTruncatedLc(r, intersection);
            return GetCorrelation(llc, rlc);
        }

        static public double GetImsCorrelation(Feature l, Feature r)
        {
            if (l == null || r == null) return -1;
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
                lc[i-sp] = feature.LcApexPeakProfile[i];
            return lc;
        }

        static private double[] GetTruncatedIms(Feature feature, Rectangle intersection)
        {
            var sp = intersection.Top - feature.ScanImsStart;
            var ep = sp + intersection.Height;
            var ims = new double[intersection.Height];
            for (var i = sp; i < ep; i++)
                ims[i-sp] = feature.ImsApexPeakProfile[i];
            return ims;
        }

        static public double GetIsotopeCorrelation(Feature[] f, double[] i)
        {
            var j = new double[f.Length];
            var max = 0.1f;
            for (var k = 0; k < j.Length; k++)
            {
                if(f[k]!=null) max = Math.Max(max,f[k].IntensityMax);
            }
            for (var k = 0; k < j.Length;k++ )
            {
                j[k] = f[k] == null? 0.1 : f[k].IntensityMax/max;
            }
            return GetCorrelation(j, i);
        }
        
        static public double GetCorrelation(double[] v1, double[] v2)
        {
            if (v1.Length <= 1) return 0.0;

           /* var d = 0.0;
            var n = new double[2];

            for (var i = 0; i < v1.Length; i++)
            {
                n[0] += v1[i]*v1[i];
                n[1] += v2[i]*v2[i];
                d += v1[i]*v2[i];
            }

            if (n[0] <= 0 || n[1] <= 0) return 0;
            return d/Math.Sqrt(n[0]*n[1]);
            //*/
            var m1 = GetSampleMean(v1);
            var m2 = GetSampleMean(v2);
            var s1 = GetSampleVariance(v1, m1);
            var s2 = GetSampleVariance(v2, m2);
            if (s1 <= 0 || s2 <= 0) return 0;
            var div = Math.Sqrt(s1*s2);
            var rho = v1.Select((t, i) => (float)((t - m1) * (v2[i] - m2) / div)).Sum();
            return Math.Min(Math.Max(0, rho / (v1.Length - 1)), 1);
            //*/
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
