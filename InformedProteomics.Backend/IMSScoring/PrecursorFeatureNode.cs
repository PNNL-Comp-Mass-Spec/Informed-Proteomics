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
            var lcScore = SubScoreFactory.GetIsotopeLCCorrelationScore(LCCorrelation, Parameter);
            var imsScore = SubScoreFactory.GetIsotopeIMSCorrelationScore(IMSCorrelation, Parameter);
            var isotopeScore = SubScoreFactory.GetIsotopeIntensityCorrelationScore(IsotopeCorrelation, Parameter);
            return lcScore + imsScore + isotopeScore;
        }
    }
}
