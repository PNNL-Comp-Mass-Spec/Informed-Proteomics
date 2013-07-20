
using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureNode : FeatureNode
    {
        private readonly SubScoreFactory _scoringParams;

        public FragmentFeatureNode(IsotopomerFeatures isotopomerFeatures, IonType fragmentIonClassBase,
            GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, fragmentIonClassBase, parameter)
        {
            _scoringParams = scoringParams;
            //if (Feature != null) Console.WriteLine(this + " " + fragmentIonClassBase.Name + " " + Feature);
            
        }

        public override sealed double GetScore()
        {
            if (IsScoreCalculated) return Score;
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(FragmentIonClassBase, IsotopeCorrelation, GroupParameter);
            //var lcScore = 0.0;
            //var imsScore = 0.0;
            //if (LcCorrelation >= 0)
            //{
            //    lcScore = _scoringParams.GetIsotopeLcCorrelationScore(FragmentIonClassBase, LcCorrelation, GroupParameter);
            //    imsScore = _scoringParams.GetIsotopeImsCorrelationScore(FragmentIonClassBase, ImsCorrelation, GroupParameter);
            //}
            //Score = lcScore + imsScore + isotopeScore;
            Score = isotopeScore;
            IsScoreCalculated = true;
            return Score;
        }
    }
}
