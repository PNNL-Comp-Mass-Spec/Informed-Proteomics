using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring.Scoring
{
    public class RankScore
    {
        public double Mass { get; set; }
        public RankScore(string fileName)
        {
            _trainingParameters = new TrainingParameters(fileName);
            Mass = 0.0;
        }

        public double GetScore(IonType ionType, int rankNum, int charge)
        {
            var targetRankTable = _trainingParameters.GetRankTable(charge, Mass, false);
            var decoyRankTable = _trainingParameters.GetRankTable(charge, Mass, true);
            var prob = targetRankTable.GetRankProbability(rankNum, ionType).Found;
            var decoyProb = decoyRankTable.GetRankProbability(rankNum, ionType).Found;

            //if (prob == 0) Console.WriteLine("*** Target: {0} {1}", ionType.Name, rankNum);
            //if (decoyProb == 0) Console.WriteLine("*** Decoy: {0} {1}", ionType.Name, rankNum);
            if (Math.Abs(prob) < 1e-10 || Math.Abs(decoyProb) < 1e-10) return 0;
            return (Math.Log(prob) - Math.Log(decoyProb));
        }

        public IonType[] GetIonTypes(int charge)
        {
            return _trainingParameters.GetIonTypes(charge, Mass);
        }

        private readonly TrainingParameters _trainingParameters;
    }
}
