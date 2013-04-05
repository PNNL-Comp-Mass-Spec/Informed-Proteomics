using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class ImsScorer
    {
        private readonly ImsDataCached _imsData;
        private readonly Ion _precursorIon;
       
        public ImsScorer(ImsDataCached imsData, Ion precursorIon) // precursorComposition has protons
        {
            _imsData = imsData;
            _precursorIon = precursorIon;
        }

        public double GetCutScore(char nTermAA, char cTermAA, Composition cutComposition, Feature precursorFeature)
        {
            var parameter = new GroupParameter(cutComposition, nTermAA, cTermAA, _precursorIon);
            return new FragmentFeatureGraph(_imsData, precursorFeature, _precursorIon, cutComposition, parameter).Score;
        }
    }
}
