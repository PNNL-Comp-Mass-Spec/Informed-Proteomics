using System;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.Util;

namespace InformedProteomics.FeatureFinding.Alignment
{
    public class LcMsFeatureAlignComparer : INodeComparer<LcMsFeature>
    {
        public LcMsFeatureAlignComparer(Tolerance tolerance, bool oneDaltonShift = false)
        {
            _oneDaltonShift = oneDaltonShift;
            _tolerance = tolerance;
        }

        private readonly bool _oneDaltonShift;
        private readonly Tolerance _tolerance;

        public bool SameCluster(LcMsFeature f1, LcMsFeature f2)
        {
            if (f1.DataSetId == f2.DataSetId)
            {
                return false;
            }

            // tolerant in mass dimension?
            if (!_oneDaltonShift)
            {
                var massTol = Math.Min(_tolerance.GetToleranceAsMz(f1.Mass), _tolerance.GetToleranceAsMz(f2.Mass));
                if (Math.Abs(f1.Mass - f2.Mass) > massTol)
                {
                    return false;
                }
            }
            else
            {
                var massTol = Math.Min(_tolerance.GetToleranceAsMz(f1.Mass), _tolerance.GetToleranceAsMz(f2.Mass));
                var massDiff = Math.Abs(f1.Mass - f2.Mass);

                if (f1.Mass > 10000 && f2.Mass > 10000)
                {
                    if (massDiff > massTol && Math.Abs(massDiff - 1) > massTol && Math.Abs(massDiff - 2) > massTol)
                    {
                        return false;
                    }
                }
                else
                {
                    if (massDiff > massTol && Math.Abs(massDiff - 1) > massTol)
                    {
                        return false;
                    }
                }
            }
            /*
            var coeLen = f1.CoElutionNetLength(f2);
            if (coeLen > f1.NetLength * 0.25 || coeLen > f2.NetLength * 0.25) return true;

            // tolerant in elution time dimension?
            var lenDiff = Math.Abs(f1.NetLength - f2.NetLength) / Math.Max(f1.NetLength, f2.NetLength);
            if (lenDiff > 0.8) return false;
            */
            //if (f1.CoElutedByNet(f2, 0.01)) return true; //e.g) 200*0.001 = 0.2 min = 30 sec

            if (f1.CoElutedByNet(f2, 0.01))
            {
                return true; //e.g) 200*0.001 = 0.2 min = 30 sec
            }

            //if (NetDiff(f1, f2) < TolNet) return true;
            return false;
        }
    }
}
