using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public GroupParameter Parameter { get; private set; }
        public Feature Feature { get; private set; }
        internal double Score { get; set; }

        protected FeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
        {
            Feature = isotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(0);
            Parameter = parameter;

            //var c = isotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(0);
            //var l = isotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(-1);
            //var r = isotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(1);
            //var r2 = isotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(2);

            //_lcCorrelation = StatisticsTools.GetLCCorrelation(c, r) - StatisticsTools.GetLCCorrelation(l, c);
            //_imsCorrelation = StatisticsTools.GetIMSCorrelation(c, r) - StatisticsTools.GetIMSCorrelation(l, c);
            
            //_isotopeCorrelation = StatisticsTools.GetCorrelation()

           
        }
        internal abstract double GetScore();
    }
}
