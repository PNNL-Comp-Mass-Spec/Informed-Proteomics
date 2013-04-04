using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FragmentFeatureNode : FeatureNode
    {
        public IonType FragmentIonClassBase { get; private set; }
        private readonly float _isotopeCorrelation, _lcCorrelation, _imsCorrelation;
        
        public FragmentFeatureNode(IsotopomerFeatures isotopomerFeatures, IonType fragmentIonClassBase, GroupParameter parameter)
            : base(isotopomerFeatures, parameter)
        {
            FragmentIonClassBase = fragmentIonClassBase;
            Score = GetScore();
        }

        internal override sealed double GetScore()
        {
            throw new System.NotImplementedException();
        }
    }
}
