using System;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureNode
    {
        public GroupParameter Parameter { get; private set; }
        public Feature Feature { get; private set; }
        internal float Score { get; set; }
        public FeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
        {
            Feature = isotopomerFeatures.GetNthFeatureFromTheMostInstenseFeature(0);
            Parameter = parameter;
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
