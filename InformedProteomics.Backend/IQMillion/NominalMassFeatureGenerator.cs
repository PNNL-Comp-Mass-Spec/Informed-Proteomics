using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.TheorFeatureGenerator;
using DeconTools.Backend.Utilities;
using DeconTools.Backend.Utilities.IsotopeDistributionCalculation;
using DeconTools.Utilities;

namespace InformedProteomics.Backend.IQMillion
{
    class NominalMassFeatureGenerator : ITheorFeatureGenerator
    {
        private static readonly IsotopicDistributionCalculator IsotopicDistCalculator = IsotopicDistributionCalculator.Instance;
        private static readonly double LowPeakCutOff = 0.005;
        public override void LoadRunRelatedInfo(ResultCollection results)
        {
            // do nothing
        }

        public override void GenerateTheorFeature(TargetBase mt)
        {
            Check.Require(mt != null, "FeatureGenerator failed. Target must not be null.");
            Check.Require(mt.MZ != null, "FeatureGenerator failed. Target must have mz.");
            Check.Require(mt.ChargeState != null, "FeatureGenerator failed. Target must have charge state.");
            IsotopicProfile iso = IsotopicDistCalculator.GetAveraginePattern(mt.MZ);
            PeakUtilities.TrimIsotopicProfile(iso, LowPeakCutOff);
            iso.ChargeState = mt.ChargeState;
            mt.IsotopicProfile = iso;
        }
    }
}
