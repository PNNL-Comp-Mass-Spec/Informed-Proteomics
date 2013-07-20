using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        private readonly SubScoreFactory _scoringParams;
        public PrecursorFeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, null, parameter)
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
            //var lcScore = 0.0;
            //var imsScore = 0.0;
            //if (LcCorrelation >= 0)
            //{
                //Console.WriteLine("precursor" + LcCorrelation + " " + ImsCorrelation);
                //lcScore = _scoringParams.GetIsotopeLcCorrelationScore(LcCorrelation, GroupParameter);
                //imsScore = _scoringParams.GetIsotopeImsCorrelationScore(ImsCorrelation, GroupParameter);
            //}
            //Score = lcScore + imsScore + isotopeScore;
            Score = isotopeScore;
            IsScoreCalculated = true;
            return Score;
        }
    }
}
