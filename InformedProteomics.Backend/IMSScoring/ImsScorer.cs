using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class ImsScorer
    {
        private readonly ImsDataCached _imsData;
        private readonly Composition _precursorComposition;
        private readonly int _charge;

        public ImsScorer(ImsDataCached imsData, Composition precursorComposition, int charge) // will precursorComposition have protons?
        {
            _imsData = imsData;
            _precursorComposition = precursorComposition;
            _charge = charge;
        }

        public double GetCutScore(char nTermAA, char cTermAA, Composition cutComposition, Feature precursorFeature)
        {
            var parameter = new GroupParameter(cutComposition, nTermAA, cTermAA, _precursorComposition, _charge);
            return new FragmentFeatureGraph(_imsData, precursorFeature, _precursorComposition, cutComposition, parameter).Score;
        }
    }
}
