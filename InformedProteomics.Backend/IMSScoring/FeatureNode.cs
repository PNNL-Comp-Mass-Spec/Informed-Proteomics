using System;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureNode
    {
        public GroupParameter Parameter { get; private set; }
        public Feature Feature { get; private set; }
        private readonly float _isotopeCorrelation, _lcCorrelation, _imsCorrelation;
        internal float Score { get; set; }
        public FeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
        {
            Feature = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(0);
            Parameter = parameter;

            var c = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(0);
            var l = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(-1);
            var r = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(1);
            var r2 = isotopomerFeatures.GetNthFeatureFromTheMostIntenseFeature(2);

            _lcCorrelation = StatisticsTools.GetLCCorrelation(c, r) - StatisticsTools.GetLCCorrelation(l, c);
            _imsCorrelation = StatisticsTools.GetIMSCorrelation(c, r) - StatisticsTools.GetIMSCorrelation(l, c);
            
            //_isotopeCorrelation = StatisticsTools.GetCorrelation()

            Score = GetScore();
        }

        

        private float GetScore()
        {
            var score = 0f;
            // unshifted corr - shifted corr =  ??
            //score += SubScoreFactory.GetIsotopeCorrelationScore(_imsCorrelation);
            //score += SubScoreFactory.GetLCCorrelationScore(_lcCorrelation);
            //score += SubScoreFactory.GetIMSCorrelationScore(_imsCorrelation);
            throw new NotImplementedException();
        }

      
    }
}
