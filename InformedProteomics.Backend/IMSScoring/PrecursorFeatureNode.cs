using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        public PrecursorFeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter)
            : base(isotopomerFeatures, parameter)
        {
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            // when calculating score, only mass index and charge values are used in parameter. When writing the parameter file, only they should be written.
            var isotopeScore = SubScoreFactory.GetIsotopeIntensityCorrelationScore(IsotopeCorrelation, GroupParameter);
            var lcScore = 0.0;
            var imsScore = 0.0;
            if (Feature != null)
            {
                lcScore = SubScoreFactory.GetIsotopeLCCorrelationScore(LCCorrelation, GroupParameter);
                imsScore = SubScoreFactory.GetIsotopeIMSCorrelationScore(IMSCorrelation, GroupParameter);
            }
            return lcScore + imsScore + isotopeScore;
        }
    }
}
