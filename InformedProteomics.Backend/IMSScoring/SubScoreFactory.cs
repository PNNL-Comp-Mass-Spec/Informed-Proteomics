using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public static class SubScoreFactory // TODO
    {
        public static void Read()
        {
             //
        }

        public static List<IonType> GetIonTypes(FragmentParameter parameter, int charge)
        {
            return null;
        }

        internal static float GetRatioScore(IonType ionType1, IonType ionType2, int ratio, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetLCCorrelationScore(IonType ionType1, IonType ionType2, float lcCorr, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetIMSCorrelationScore(IonType ionType1, IonType ionType2, float imsCorr, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        public static float GetKLDivergence(IonType ionType, IonType ionType1, int ratio, float lcCorrelation, float imsCorrelation, FragmentParameter parameter)
        {
            throw new System.NotImplementedException();
        }
    }
}
