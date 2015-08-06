using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassFeature
{
    public interface ILcMsFeatureComparer
    {
        bool Match(LcMsFeature f1, LcMsFeature f2);
    }

    public class LcMsFeatureDefaultComparer : ILcMsFeatureComparer
    {
        public LcMsFeatureDefaultComparer(Tolerance tolerance, bool oneDaltonShift = false)
        {
            _oneDaltonShift = oneDaltonShift;
            _tolerance = tolerance;
        }

        private readonly bool _oneDaltonShift;
        private readonly Tolerance _tolerance;

        public bool Match(LcMsFeature f1, LcMsFeature f2)
        {
            if (f1.DataSetId == f2.DataSetId) return false;
            // tolerant in mass dimension?
            if (!_oneDaltonShift)
            {
                var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));
                if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;
            }
            else
            {
                var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));
                var massDiff = Math.Abs(f1.Mass - f2.Mass);

                if (f1.Mass > 10000 && f2.Mass > 10000)
                {
                    if (massDiff > massTol && Math.Abs(massDiff - 1) > massTol && Math.Abs(massDiff - 2) > massTol) return false;
                }
                else
                {
                    if (massDiff > massTol && Math.Abs(massDiff - 1) > massTol) return false;
                }
            }

            // tolerant in elution time dimension?
            var lenDiff = Math.Abs(f1.NetLength - f2.NetLength) / Math.Min(f1.NetLength, f2.NetLength);
            if (lenDiff > 0.5) return false;

            if (f1.CoElutedByNet(f2, 0.001)) return true; //e.g) 200*0.001 = 0.2 min = 30 sec
            //if (NetDiff(f1, f2) < TolNet) return true;
            return false;
        }
    }

    public class LcMsFeatureComparer : ILcMsFeatureComparer
    {
        public bool Match(LcMsFeature f1, LcMsFeature f2)
        {
            throw new NotImplementedException();
        }
    }
}
