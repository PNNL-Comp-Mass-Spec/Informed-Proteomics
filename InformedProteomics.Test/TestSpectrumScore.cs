using System;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.LikelihoodScoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSpectrumScore
    {
        [Test]
        public void RankScore()
        {
            var ranks = 151;
            var rankScorer = new RankScore(@"\\protoapps\UserData\Wilkins\MSGFPlusTrainingData\HCD_QExactive_Tryp\HCD_QExactive_Tryp_RankProbabilities.txt");
            for (int charge = 1; charge < 4; charge++)
            {
                var ionTypes = rankScorer.GetIonTypes(charge);
                foreach (var ionType in ionTypes)
                {
                    for (int r = 1; r < ranks; r++)
                    {
                        Console.WriteLine("Charge: {0}, Ion Type: {1}, Rank: {2}, Score: {3}",
                            charge, ionType.Name, r+1, rankScorer.GetScore(ionType, r, charge));
                    }
                }
            }
        }
    }
}
