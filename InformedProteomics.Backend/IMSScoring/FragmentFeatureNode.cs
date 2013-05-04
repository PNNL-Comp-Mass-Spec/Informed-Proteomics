
using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureNode : FeatureNode
    {
        public IonType FragmentIonClassBase { get; private set; }
        
        public FragmentFeatureNode(IsotopomerFeatures isotopomerFeatures, IonType fragmentIonClassBase, GroupParameter parameter)
            : base(isotopomerFeatures, parameter)
        {
            FragmentIonClassBase = fragmentIonClassBase;
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            var lcScore = SubScoreFactory.GetIsotopeLCCorrelationScore(FragmentIonClassBase, LCCorrelation, GroupParameter);
            var imsScore = SubScoreFactory.GetIsotopeIMSCorrelationScore(FragmentIonClassBase, IMSCorrelation, GroupParameter);
            var isotopeScore = SubScoreFactory.GetIsotopeIntensityCorrelationScore(FragmentIonClassBase, IsotopeCorrelation, GroupParameter);
            var score = lcScore + imsScore + isotopeScore;
            return score;
        }
    }
}
