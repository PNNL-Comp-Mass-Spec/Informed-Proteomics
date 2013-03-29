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

        public static List<IonType> GetIonTypes(GroupParameter parameter, int charge)
        {
            return null;
        }

        internal static float GetRatioScore(IonType ionType1, IonType ionType2, int ratio, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetLCCorrelationScore(IonType ionType1, IonType ionType2, float lcCorr, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetIMSCorrelationScore(IonType ionType1, IonType ionType2, float imsCorr, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        public static float GetKLDivergence(IonType ionType, IonType ionType1, int ratio, float lcCorrelation, float imsCorrelation, GroupParameter parameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetKLDivergence(IonType ionType, int _ratio, float _lcCorrelation, float _imsCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetRatioScore(IonType ionType, int _ratio, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetLCCorrelationScore(IonType ionType, float _lcCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetIMSCorrelationScore(IonType ionType, float _imsCorrelation, GroupParameter groupParameter)
        {
            throw new System.NotImplementedException();
        }
    }
}
