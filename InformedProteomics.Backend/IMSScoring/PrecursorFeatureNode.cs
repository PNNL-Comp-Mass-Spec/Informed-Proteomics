using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        private readonly SubScoreFactory _scoringParams;
        public PrecursorFeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, null, isotopomerFeatures.GetMostAbundantFeature(), parameter)
        {
            _scoringParams = scoringParams;
            Score = GetScore();
        }

        public override sealed double GetScore()
        {
            // when calculating score, only mass index and charge values are used in parameter. When writing the parameter file, only they should be written.
            if (IsScoreCalculated) return Score;
           // Console.WriteLine(this + " " + IsotopeCorrelation);
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(IsotopeCorrelation, GroupParameter);
            
            var lcScore = _scoringParams.GetLcCorrelationScore(LcCorrelation, GroupParameter);
            var imsScore = _scoringParams.GetImsCorrelationScore(ImsCorrelation, GroupParameter);
           
            Score = lcScore + imsScore + isotopeScore;
            //Score = isotopeScore;
            IsScoreCalculated = true;
            //Console.WriteLine(this + " " + lcScore + " " + imsScore + " " + isotopeScore);
            return Score;
        }
    }
}
