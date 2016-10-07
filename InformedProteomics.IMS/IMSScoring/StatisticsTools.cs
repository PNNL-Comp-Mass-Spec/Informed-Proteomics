using System;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.Utils;
using InformedProteomics.IMS.IMS;

namespace InformedProteomics.IMS.IMSScoring
{
    public class StatisticsTools
    {
        static public double GetLcCorrelation(Feature l, Feature r)
        {
            if (l == null || r == null) return -1;
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var llc = GetTruncatedLc(l, intersection);
            var rlc = GetTruncatedLc(r, intersection);
            return SimpleMath.GetCorrelation(llc, rlc);
        }

        static public double GetImsCorrelation(Feature l, Feature r)
        {
            if (l == null || r == null) return -1;
            var intersection = Rectangle.Intersect(l.GetBoundary(), r.GetBoundary());
            var lims = GetTruncatedIms(l, intersection);
            var rims = GetTruncatedIms(r, intersection);
            return SimpleMath.GetCorrelation(lims, rims);
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
                if(f[k]!=null) max = Math.Max(max,(float) f[k].IntensityMax);
            }
            for (var k = 0; k < j.Length;k++ )
            {
                j[k] = f[k] == null? 0.1 : f[k].IntensityMax/max;
            }
            return SimpleMath.GetCorrelation(j, i);
        }
    }
}
