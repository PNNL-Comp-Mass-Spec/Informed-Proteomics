using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureScore
    {
        public double EnvelopeCorrelation;
        public double EnvelopeCorrelationSummed;

        public double RankSum;
        public double Poisson;

        public double BhattacharyyaDistance;
        public double BhattacharyyaDistanceSummed;
        public double BhattacharyyaDistanceSummedOverCharges;
        public double BhattacharyyaDistanceSummedOverTimes;

        public double XicCorrMean;
        public double XicCorrMin;

        public double BhattacharyyaDistanceSummedOverEvenCharges;
        public double BhattacharyyaDistanceSummedOverOddCharges;

        public double AbundanceChangesOverCharges;

        public double Likelihood;
    }
}
