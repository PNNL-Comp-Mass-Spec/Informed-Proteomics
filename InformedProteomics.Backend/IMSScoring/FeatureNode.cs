using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureNode
    {
        public FragmentParameter Parameter { get; private set; }
        public IonType IonType { get; private set; }
        public Feature Feature { get; private set; }

        public FeatureNode(Feature feature, IonType ionType, FragmentParameter parameter)
        {
            Feature = feature;
            IonType = ionType;
            Parameter = parameter;
        }
    }
}
