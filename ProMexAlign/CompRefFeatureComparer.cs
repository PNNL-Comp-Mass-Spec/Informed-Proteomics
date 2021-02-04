using System;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.Util;

namespace PromexAlign
{
    public class CompRefFeatureComparer : INodeComparer<LcMsFeature>
    {
        public CompRefFeatureComparer(Tolerance tolerance = null)
        {
            _tolerance = tolerance ?? new Tolerance(10);
        }

        public bool SameCluster(LcMsFeature f1, LcMsFeature f2)
        {
            if (f1.DataSetId == f2.DataSetId)
            {
                return false;
            }
            // tolerant in mass dimension?

            var massTol = Math.Min(_tolerance.GetToleranceAsMz(f1.Mass), _tolerance.GetToleranceAsMz(f2.Mass));
            if (Math.Abs(f1.Mass - f2.Mass) > massTol)
            {
                return false;
            }

            //if (!f1.CoElutedByNet(f2, 0.004)) return false; //e.g) 200*0.001 = 0.2 min = 30 sec
            if (f1.ProteinSpectrumMatches.Count > 0 && f2.ProteinSpectrumMatches.Count > 0)
            {
                return f1.CoElutedByNet(f2, 0.03) && f1.ProteinSpectrumMatches.ShareProteinId(f2.ProteinSpectrumMatches);
            }

            return f1.CoElutedByNet(f2, 0.01);
        }

        private readonly Tolerance _tolerance;
    }
}
