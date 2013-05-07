
using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureNode : FeatureNode
    {
        public IonType FragmentIonClassBase { get; private set; }

        private readonly SubScoreFactory _scoringParams;

        public FragmentFeatureNode(IsotopomerFeatures isotopomerFeatures, IonType fragmentIonClassBase,
            GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, parameter)
        {
            FragmentIonClassBase = fragmentIonClassBase;
            _scoringParams = scoringParams;
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(FragmentIonClassBase, IsotopeCorrelation, GroupParameter);
            var lcScore = 0.0;
            var imsScore = 0.0;
            if (Feature != null)
            {
                lcScore = _scoringParams.GetIsotopeLcCorrelationScore(FragmentIonClassBase, LcCorrelation, GroupParameter);
                imsScore = _scoringParams.GetIsotopeImsCorrelationScore(FragmentIonClassBase, ImsCorrelation, GroupParameter);
            }
            var score = lcScore + imsScore + isotopeScore;
            return score;
        }
    }
}
