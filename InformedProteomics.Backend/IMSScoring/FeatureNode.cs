using System;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public GroupParameter Parameter { get; private set; }
        public Feature Feature { get; private set; }
        internal float Score { get; set; }

        protected FeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
        {
            Feature = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(0);
            Parameter = parameter;

            //var c = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(0);
            //var l = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(-1);
            //var r = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(1);
            //var r2 = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(2);

            //_lcCorrelation = StatisticsTools.GetLCCorrelation(c, r) - StatisticsTools.GetLCCorrelation(l, c);
            //_imsCorrelation = StatisticsTools.GetIMSCorrelation(c, r) - StatisticsTools.GetIMSCorrelation(l, c);
            
            //_isotopeCorrelation = StatisticsTools.GetCorrelation()

           // Score = GetScore();
        }

        internal abstract float GetScore();
    }
}
