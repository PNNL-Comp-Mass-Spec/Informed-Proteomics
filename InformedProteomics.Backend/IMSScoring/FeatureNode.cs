using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureNode
    {
        public FragmentParameter Parameter { get; private set; }
        public IonType FragmentIonClassBase { get; private set; }
        public Feature Feature { get; private set; }

        public FeatureNode(Feature feature, IonType fragmentIonClassBase, FragmentParameter parameter)
        {
            Feature = feature;
            FragmentIonClassBase = fragmentIonClassBase;
            Parameter = parameter;
        }
    }
}
