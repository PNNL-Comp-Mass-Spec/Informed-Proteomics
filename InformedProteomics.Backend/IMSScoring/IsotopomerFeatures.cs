using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class IsotopomerFeatures : List<Feature>
    {
        public IsotopomerFeatures(ImsDataCached imsData, Composition composition)
        {
            //TODO
        }

        public Feature GetNthFeatureFromTheMostInstenseFeature(int n)
        {
            throw new NotImplementedException();
        }
    }
}
