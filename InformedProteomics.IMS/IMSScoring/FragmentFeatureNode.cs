using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.IMS.IMS;

namespace InformedProteomics.IMS.IMSScoring
{
    public class FragmentFeatureNode : FeatureNode
    {
        private readonly SubScoreFactory _scoringParams;

        public FragmentFeatureNode(IsotopomerFeatures isotopomerFeatures, IonType fragmentIonClassBase, Feature precursorFeature,
            GroupParameter parameter, SubScoreFactory scoringParams)
            : base(isotopomerFeatures, fragmentIonClassBase, precursorFeature, parameter)
        {
            _scoringParams = scoringParams;
        }

        public override sealed double GetScore()
        {
            if (IsScoreCalculated) return Score;
            var isotopeScore = _scoringParams.GetIsotopeIntensityCorrelationScore(IsotopeCorrelation, GroupParameter);

            //var lcScore = _scoringParams.GetLcCorrelationScore(LcCorrelation, GroupParameter);
            //var imsScore = _scoringParams.GetImsCorrelationScore(ImsCorrelation, GroupParameter);

            //Score = lcScore + imsScore + isotopeScore;
            Score = isotopeScore;
            //Console.WriteLine(this + " " + lcScore + " " + imsScore + " " + isotopeScore);

            IsScoreCalculated = true;
            return Score;
        }
    }
}
