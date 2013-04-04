using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class PrecursorFeatureNode : FeatureNode
    {
        public PrecursorFeatureNode(Feature precursorFeature, GroupParameter parameter)
            : base(precursorFeature, parameter)
        {
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            return 0;
        }
    }
}
