using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class IsotopomerFeatures : List<Feature>
    {
        public IsotopomerFeatures(ImsDataCached imsData, Composition composition, Feature precursorFeature, int charge)
        {
                       
            //TODO
        }

        public Feature GetNthFeatureFromTheMostIntenseFeature(int n)
        {
            throw new NotImplementedException();
        }

        public float GetIntensityRatioOfNthFeatureFromTheMostIntenseFeature(int n)
        {
            throw new NotImplementedException();            
        }
    }
}
