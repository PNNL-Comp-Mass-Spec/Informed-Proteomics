using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public const int NumSupport = 4; // used for training
        public const int NumMinusIsotope = 1; // used for training
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
            if (Feature.IntensityMax <= 0) return;
            var f = new Feature[NumSupport];
            var i = new double[NumSupport];
            for (var k = 0; k < NumSupport; k++)
            {
                f[k] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(k - NumMinusIsotope);
                i[k] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(k - NumMinusIsotope);
            }
            LCCorrelation = StatisticsTools.GetLCCorrelation(f[NumMinusIsotope], f[NumMinusIsotope + 1]);
            IMSCorrelation = StatisticsTools.GetIMSCorrelation(f[NumMinusIsotope], f[NumMinusIsotope + 1]);
            IsotopeCorrelation = StatisticsTools.GetIsotopeCorrelation(f, i);
        }

        internal abstract double GetScore();
    }
}
