using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public static class ScoringParameter // TODO
    {
        public static void Read()
        {
             //
        }

        public static List<IonType> GetIonTypes(FragmentParameter parameter, int charge)
        {
            return null;
        }

        internal static float GetRatioScore(IonType fragmentIonClassBaseType1, IonType fragmentIonClassBaseType2, int ratio, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetLCCorrelationScore(IonType fragmentIonClassBaseType1, IonType fragmentIonClassBaseType2, float lcCorr, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        internal static float GetIMSCorrelationScore(IonType fragmentIonClassBaseType1, IonType fragmentIonClassBaseType2, float imsCorr, FragmentParameter fragmentParameter)
        {
            throw new System.NotImplementedException();
        }

        public static float GetKLDivergence(IonType fragmentIonClassBase, IonType fragmentIonClassBaseType1, int ratio, float lcCorrelation, float imsCorrelation, FragmentParameter parameter)
        {
            throw new System.NotImplementedException();
        }
    }
}
