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
        private Feature _previousPrecursorFeature; // only for speed-up
        public PrecursorFeatureNode PrecursorFeatureNode { get; private set; }

        static public void ReadScoringParameterFile(string fileName)
        {
            SubScoreFactory.Read(fileName);
        }

        public ImsScorer(ImsDataCached imsData, Ion precursorIon) // precursorComposition has protons
        {
            _imsData = imsData;
            _precursorIon = precursorIon;
        }

        public double GetCutScore(char nTermAA, char cTermAA, Composition cutComposition, Feature precursorFeature)
        {
            var parameter = new GroupParameter(cutComposition, nTermAA, cTermAA, _precursorIon);
            if (_previousPrecursorFeature == null || _previousPrecursorFeature != precursorFeature)
            {
                PrecursorFeatureNode =
                    new PrecursorFeatureNode(
                        IsotopomerFeatures.GetPrecursorIsotopomerFeatures(_imsData, _precursorIon, precursorFeature),
                        parameter);
                _previousPrecursorFeature = precursorFeature;
            }
            return new FragmentFeatureGraph(_imsData, PrecursorFeatureNode, precursorFeature, _precursorIon, cutComposition, parameter).Score;
        }
    }
}
