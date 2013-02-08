using System;

namespace InformedProteomics.Backend.Scoring
{
    public static class SubScoreFactory
    {
        public static void Read(string fileName)
        {
        }

        private static float GetNormalLogRatio(float x, float numVar)
        {
            x += (numVar-1)*.3f;
            var sigma = (float)Math.Sqrt(0.5);
            const float mean1 = 0.9f;
            const float mean2 = 0.2f;
            return (float)(-.5f * (Math.Pow((x - mean1) / sigma, 2) - Math.Pow((x - mean2) / sigma, 2)));
        }

        internal static float GetPrecursorIonLikelihoodRatioScore(float rawScore, float numVar) // spectrum para
        {
            return GetNormalLogRatio(rawScore, numVar);
        }

        internal static float GetProductIonXICLikelihoodRatioScore(float rawScore, float numVar) // spectrum para
        {
            return GetNormalLogRatio(rawScore, numVar);
        }

        internal static float GetPrecursorIonCorrelationCoefficient(int c1, int c2) // spectrum para
        {
            if (c1 == c2) return 1;
            return 0.8f;
        }

        internal static float GetProductIonCorrelationCoefficient(IonType ion1, IonType ion2) // spectrum para, peak para (of ion1)
        {
            if (ion1.Equals(ion2)) return 1;
            if (ion1.Charge != ion2.Charge) return 0.4f;
            if (ion1.IsPrefix == ion2.IsPrefix) return 0.7f;
           
            return 0.5f;
            
        }

        internal static float GetProductIonCorrelationCoefficient(IonType ion) // spectrum para, peak para
        {
            return 0.5f;
        }

        public static float GetProductIonSpectrumScore(IonType ion1, IonType ion2, double intensity1, double intensity2)
        {
            float rawScore;
            if (ion1 == null)
            {
                if (ion2.Type.Equals("y")) rawScore = 0.7f / 0.1f;
                else if (ion2.Type.Equals("b")) rawScore = 0.6f / 0.1f;
                else if (ion2.Type.Equals("a")) rawScore = 0.3f/0.1f;
                else rawScore = 0.2f/0.1f;
            }else
            {
                rawScore = GetRatioScore(intensity1, intensity2);
            }
            return (float) Math.Log(rawScore);
        }

        private static float GetRatioScore(double i1, double i2)
        {
            const float noiseProb = 0.1f/5.0f;
            float rawScore;

            if (i1 <= 0)
            {
                rawScore = 0.01f/noiseProb;
            }
            else
            {
                var r = (float)(i2/i1);
                if (r > 2) rawScore = 0.04f/noiseProb;
                else if (r > 1) rawScore = 0.05f/noiseProb;
                else if (r > .5) rawScore = 0.6f/noiseProb;
                else if (r > 0) rawScore = 0.2f/noiseProb;
                else rawScore = 0.1f/(1-noiseProb*5);

            }

            return rawScore;
        }

    }
}
