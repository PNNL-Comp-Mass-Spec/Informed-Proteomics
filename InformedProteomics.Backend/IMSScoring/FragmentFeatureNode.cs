
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
            
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(FragmentIonClassBase, IsotopeCorrelation, GroupParameter);
            var lcScore = 0.0;
            var imsScore = 0.0;
            if (LcCorrelation >= 0)
            {
                lcScore = _scoringParams.GetIsotopeLcCorrelationScore(FragmentIonClassBase, LcCorrelation, GroupParameter);
                imsScore = _scoringParams.GetIsotopeImsCorrelationScore(FragmentIonClassBase, ImsCorrelation, GroupParameter);
            }
            var score = lcScore + imsScore + isotopeScore;
            //if (Feature != null) Console.WriteLine("Ion Type " + FragmentIonClassBase.Name);
            
            return score;
        }
    }
}
