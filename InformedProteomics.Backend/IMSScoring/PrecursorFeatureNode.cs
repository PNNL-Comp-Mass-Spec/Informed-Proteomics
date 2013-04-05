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
            return 0;
        }
    }
}
