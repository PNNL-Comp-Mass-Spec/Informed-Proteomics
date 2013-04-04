using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public static class SubScoreFactory // TODO
    {
        public static void Read()
        {
             //
        }

        public static List<IonType> GetIonTypes(GroupParameter parameter)
        {
            return null;
        }

        internal static double GetRatioScore(IonType ionType1, IonType ionType2, int ratio, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetLCCorrelationScore(IonType ionType1, IonType ionType2, double lcCorr, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetIMSCorrelationScore(IonType ionType1, IonType ionType2, double imsCorr, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        public static double GetKLDivergence(IonType ionType, IonType ionType1, int ratio, double lcCorrelation, double imsCorrelation, GroupParameter parameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetKLDivergence(IonType ionType, int _ratio, double _lcCorrelation, double _imsCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetRatioScore(IonType ionType, int _ratio, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetLCCorrelationScore(IonType ionType, double _lcCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetIMSCorrelationScore(IonType ionType, double _imsCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }
    }
}
