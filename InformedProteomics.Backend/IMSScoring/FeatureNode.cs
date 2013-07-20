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
        public IsotopomerFeatures IsotopomerFeatures { get; private set; }
        public Feature Feature { get; private set; }
        public IonType FragmentIonClassBase { get; private set; }
        protected double Score { get; set; }
        protected bool IsScoreCalculated = false; // for speed-up
        internal double IsotopeCorrelation;//, LcCorrelation, ImsCorrelation;
        
        protected FeatureNode(IsotopomerFeatures isotopomerFeatures, IonType ionType, GroupParameter groupParameter)
        {
            IsotopomerFeatures = isotopomerFeatures;
            FragmentIonClassBase = ionType;
            //Feature = IsotopomerFeatures.GetMostAbundantFeature();
            Feature = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(0);
            GroupParameter = groupParameter;
            GetCorrelations();
        }

        private void GetCorrelations()
        {
            var f = new Feature[IsotopomerFeatures.Count];
            var i = new double[IsotopomerFeatures.Count];
            //Console.WriteLine(this + " " + f.Length);
            for (var k = 0; k < f.Length; k++)
            {
                f[k] = IsotopomerFeatures[k];
                i[k] = IsotopomerFeatures.TheoreticalIsotopomerEnvelope[k];
             //   Console.WriteLine(k + " " + (f[k] == null ? 0 : f[k].IntensityMax) + " " + i[k]);
            }
          //  LcCorrelation = StatisticsTools.GetLcCorrelation(f[OffsetFromMonoIsotope], f[OffsetFromMonoIsotope + 1]);
          //  ImsCorrelation = StatisticsTools.GetImsCorrelation(f[OffsetFromMonoIsotope], f[OffsetFromMonoIsotope + 1]);
            IsotopeCorrelation = StatisticsTools.GetIsotopeCorrelation(f, i);
            //Console.WriteLine(this + " Node " + LcCorrelation + " "+ ImsCorrelation + " " + IsotopeCorrelation);
        }

        public abstract double GetScore();
    }
}
