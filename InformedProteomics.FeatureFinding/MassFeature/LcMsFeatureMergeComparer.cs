using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public class LcMsFeatureMergeComparer : INodeComparer<LcMsPeakCluster>
    {
        public LcMsFeatureMergeComparer(Tolerance tolerance)
        {
            _tolerance = tolerance;
        }
            
        public bool SameCluster(LcMsPeakCluster f1, LcMsPeakCluster f2)
        {
            var massTh = _tolerance.GetToleranceAsTh(f1.RepresentativeMass);
            var massDiff = Math.Abs(f1.RepresentativeMass - f2.RepresentativeMass);
            if (massDiff > massTh) return false;
            
            // close in elution time
            if (f1.CoElutedByNet(f2, 0.01)) return true;
            return false;
        }
        private readonly Tolerance _tolerance;
    }
    /*
    public class LcMsProteoformMergeCompare : INodeComparer<LcMsPeakCluster>
    {
        public LcMsProteoformMergeCompare(Tolerance tolerance)
        {
            _tolerance = tolerance;
        }

        public bool SameCluster(LcMsPeakCluster f1, LcMsPeakCluster f2)
        {
            var massTh = _tolerance.GetToleranceAsTh(f1.RepresentativeMass);
            var massDiff = Math.Abs(f1.RepresentativeMass - f2.RepresentativeMass);
            if (massDiff > massTh) return false;

            // significantly overlapped features
            var coeLen = f1.CoElutionLength(f2);
            if (coeLen > f2.ElutionLength * 0.6 || coeLen > f1.ElutionLength * 0.6) return true;

            // close in elution time
            if (f1.CoElutedByNet(f2, 0.003))
            {
                // two differnt hills
                if (Math.Abs(f1.ApexElutionTime - f2.ApexElutionTime) / f1.Run.GetElutionTime(f1.Run.MaxLcScan) > 0.01d &&
                    f1.ApexIntensity / f1.BoundaryIntensity > 3.0d &&
                    f2.ApexIntensity / f2.BoundaryIntensity > 3.0d)
                {
                    //Console.WriteLine("");
                    //Console.WriteLine("{0:0.00}\t{1:0.00}\t{2:0.00}\t{3:0.00}\t{4:0.00}\t{5:0.00}", f1.MinElutionTime, f1.MaxElutionTime, f1.ApexElutionTime, f1.ApexIntensity, f1.BoundaryIntensity, f1.Abundance);
                    //Console.WriteLine("{0:0.00}\t{1:0.00}\t{2:0.00}\t{3:0.00}\t{4:0.00}\t{5:0.00}", f2.MinElutionTime, f2.MaxElutionTime, f2.ApexElutionTime, f2.ApexIntensity, f2.BoundaryIntensity, f2.Abundance);
                    //Console.WriteLine("");
                    return false;
                }

                // otherwise, they are fragmentized features, which
                return true;
            }

            return false;
        }
        private readonly Tolerance _tolerance;
    }*/

}
