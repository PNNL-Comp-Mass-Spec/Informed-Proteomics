using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public GroupParameter Parameter { get; private set; }
        public IsotopomerFeatures IsotopomerFeatures { get; private set; }
        public Feature Feature { get; private set; }
        internal double Score { get; set; }
        internal double IsotopeCorrelation, LCCorrelation, IMSCorrelation;
        
        protected FeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
        {
            IsotopomerFeatures = isotopomerFeatures;
            Feature = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(0);
            Parameter = parameter;
            GetCorrelations();
        }

        private void GetCorrelations()
        {
            var f = new Feature[4];

            f[0] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(-1);
            f[1] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(0);
            f[2] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(1);
            f[3] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(2);

            var i = new double[4];
            i[0] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(-1);
            i[1] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(0);
            i[2] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(1);
            i[3] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(2);

            LCCorrelation = StatisticsTools.GetLCCorrelation(f[1], f[2]);
            IMSCorrelation = StatisticsTools.GetIMSCorrelation(f[1], f[2]);
            IsotopeCorrelation = StatisticsTools.GetIsotopeCorrelation(f, i);
        }

        internal abstract double GetScore();
    }
}
