
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
            var lcScore = SubScoreFactory.GetIsotopeLCCorrelationScore(FragmentIonClassBase, LCCorrelation, Parameter);
            var imsScore = SubScoreFactory.GetIsotopeIMSCorrelationScore(FragmentIonClassBase, IMSCorrelation, Parameter);
            var isotopeScore = SubScoreFactory.GetIsotopeIntensityCorrelationScore(FragmentIonClassBase, IsotopeCorrelation, Parameter);
            return lcScore + imsScore + isotopeScore;
        }
    }
}
