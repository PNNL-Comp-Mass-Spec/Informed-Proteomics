using System;
using InformedProteomics.FeatureFinding.Util;

namespace InformedProteomics.FeatureFinding.IsotopicEnvelope
{
    public class LcMsEnvelopeComparer : INodeComparer<ObservedIsotopeEnvelope>
    {
        public bool SameCluster(ObservedIsotopeEnvelope e1, ObservedIsotopeEnvelope e2)
        {
            //var massDiff = Math.Abs(e1.MonoMass - e2.MonoMass);
            //var tolerance = new Tolerance(5);
            //if (massDiff > tolerance.GetToleranceAsTh(e1.MonoMass)) return false;

            //var summedEnvelop = new double[e1.Size];
            //e1.SumEnvelopeTo(summedEnvelop);
            //e2.SumEnvelopeTo(summedEnvelop);

            var scanDiff = Math.Abs(e2.ScanNum - e1.ScanNum);
            if (scanDiff > 3) return false;

            var chargeDiff = Math.Abs(e2.Charge - e1.Charge);
            if (chargeDiff > 3) return false;

            return true;
        }


    }
}
