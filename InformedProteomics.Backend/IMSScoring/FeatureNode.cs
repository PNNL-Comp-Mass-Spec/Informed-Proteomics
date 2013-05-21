using System;
using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public abstract class FeatureNode
    {
        public const int NumMaxSupport = 7; // used for training
        public const int NumMinusIsotope = 1; // used for training
        public GroupParameter GroupParameter { get; private set; }
        public IsotopomerFeatures IsotopomerFeatures { get; private set; }
        public Feature Feature { get; private set; }
        public IonType FragmentIonClassBase { get; private set; }
        internal double Score { get; set; }
        internal double IsotopeCorrelation, LcCorrelation, ImsCorrelation;
        
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
            var f = new Feature[IsotopomerFeatures.SupportSize];
            var i = new double[IsotopomerFeatures.SupportSize];
            //Console.WriteLine(this + " " + f.Length);
            for (var k = 0; k < f.Length; k++)
            {
                f[k] = IsotopomerFeatures.GetNthFeatureFromTheoreticallyMostIntenseFeature(k - NumMinusIsotope);
                i[k] = IsotopomerFeatures.GetTheoreticalIntensityOfNthFeature(k - NumMinusIsotope);
            }
            LcCorrelation = StatisticsTools.GetLcCorrelation(f[NumMinusIsotope], f[NumMinusIsotope + 1]);
            ImsCorrelation = StatisticsTools.GetImsCorrelation(f[NumMinusIsotope], f[NumMinusIsotope + 1]);
            IsotopeCorrelation = StatisticsTools.GetIsotopeCorrelation(f, i);
           // Console.WriteLine(this + " Node " + LcCorrelation + " "+ ImsCorrelation + " " + IsotopeCorrelation);
        }

        internal abstract double GetScore();
    }
}
