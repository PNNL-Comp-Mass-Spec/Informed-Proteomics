using System.Collections.Generic;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        public int Charge { get;private set; }
        public PrecursorFeatureNode(IsotopomerFeatures isotopomerFeatures, GroupParameter parameter) : base(isotopomerFeatures, parameter)
        {
            Charge = parameter.Charge;
        }
    }
}
