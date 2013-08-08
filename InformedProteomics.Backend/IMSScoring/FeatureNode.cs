using System;
using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public const int IsotopeWindowSizeInDa = 6; // used for training
        public const int OffsetFromMonoIsotope = 1; // used for training
        public GroupParameter GroupParameter { get; private set; }
        private readonly IsotopomerFeatures _isotopomerFeatures;
        public Feature Feature { get; private set; }
        public Feature PrecursorFeature { get; private set; }
        public IonType FragmentIonClassBase { get; private set; }
        protected double Score { get; set; }
        protected bool IsScoreCalculated = false; // for speed-up
        internal double IsotopeCorrelation;
        public double LcCorrelation { get; private set; }
        public double ImsCorrelation { get; private set; }

        protected FeatureNode(IsotopomerFeatures isotopomerFeatures, IonType ionType, Feature precursorFeature, GroupParameter groupParameter)
        {
            _isotopomerFeatures = isotopomerFeatures;
            FragmentIonClassBase = ionType;
            PrecursorFeature = precursorFeature;
            GroupParameter = groupParameter;
            GetCorrelations();
        }

        private void GetCorrelations()
        {
            var i = new double[_isotopomerFeatures.Count];
            var f = new Feature[_isotopomerFeatures.Count];
            LcCorrelation = ImsCorrelation = .0;
            for (var k = 0; k < _isotopomerFeatures.Count; k++)
            {
                i[k] = _isotopomerFeatures.TheoreticalIsotopomerEnvelope[k];
                
                var lcCorrelation = StatisticsTools.GetLcCorrelation(PrecursorFeature, _isotopomerFeatures[k]);
                var imsCorrelation = StatisticsTools.GetImsCorrelation(PrecursorFeature, _isotopomerFeatures[k]);
                
                if (lcCorrelation >= 0.25 && imsCorrelation >= 0.25) //TODO take the good numbers!!
                {
                    f[k] = _isotopomerFeatures[k];
                }
                if (k != _isotopomerFeatures.MaxIntensityIndex) continue;
                LcCorrelation = lcCorrelation;
                ImsCorrelation = imsCorrelation;
               // Console.WriteLine((FragmentIonClassBase == null? "p" : FragmentIonClassBase.Name) + "** " + k + " " + (f[k] == null));
                Feature = f[k]; // if lc correlation score and ims correlation score are too low, Feature = null;
            }
            IsotopeCorrelation = StatisticsTools.GetIsotopeCorrelation(f, i);
        }

        public abstract double GetScore();
    }
}
