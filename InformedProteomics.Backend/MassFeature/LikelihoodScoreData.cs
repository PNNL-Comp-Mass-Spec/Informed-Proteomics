namespace InformedProteomics.Backend.MassFeature
{
    internal static class LikelihoodScoreData
    {
        internal readonly static double[][] AbuScore;
        internal readonly static double[][] CorrScore;
        internal readonly static double[][] DistScore;
        internal readonly static double[][] IntScore;

        internal readonly static double[][] SummedCorrScore;
        internal readonly static double[][] SummedDistScore;
        internal readonly static double[][] SummedIntScore;

        internal readonly static double[][] XicCorrScore1;
        internal readonly static double[][] XicCorrScore2;
        private const int Nrow = 28;
        private const int Ncol = 1001;

        static LikelihoodScoreData()
        {

           
        }
    }
}
