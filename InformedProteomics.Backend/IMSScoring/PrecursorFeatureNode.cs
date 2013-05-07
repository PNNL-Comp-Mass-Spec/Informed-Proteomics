using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        private readonly SubScoreFactory _scoringParams;
        public PrecursorFeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, parameter)
        {
            _scoringParams = scoringParams;
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            // when calculating score, only mass index and charge values are used in parameter. When writing the parameter file, only they should be written.
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(IsotopeCorrelation, GroupParameter);
            var lcScore = 0.0;
            var imsScore = 0.0;
            if (Feature != null)
            {
                lcScore = _scoringParams.GetIsotopeLcCorrelationScore(LcCorrelation, GroupParameter);
                imsScore = _scoringParams.GetIsotopeImsCorrelationScore(ImsCorrelation, GroupParameter);
            }
            return lcScore + imsScore + isotopeScore;
        }
    }
}
