using System;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring.Scoring
{
    public class RankScore
    {
        public RankScore(ActivationMethod activationMethod, Enzyme enzyme)
        {
            if (activationMethod == ActivationMethod.HCD && enzyme == Enzyme.Trypsin)
            {
                var paramFile = Properties.Resources.HCD_Trypsin;
                var stream = new MemoryStream();
                var writer = new StreamWriter(stream);
                writer.Write(paramFile);
                writer.Flush();
                stream.Position = 0;
                _trainingParameters = new TrainingParameters(stream);
            }
            else
            {
                throw new ArgumentException("No parameter file available for selected arguments.");
            }
        }

        [Obsolete("Use the constructor with two parameters: activationMethod and enzyme")]
        public RankScore(ActivationMethod activationMethod, Ms2DetectorType ms2DetectorType, Enzyme enzyme, Protocol protocol) :
            this(activationMethod, enzyme)
        {
        }

        public RankScore(string fileName)
        {
            _trainingParameters = new TrainingParameters(fileName);
        }

        public double GetScore(IonType ionType, int rankNum, int charge, double mass)
        {
            var targetRankTable = _trainingParameters.GetRankTable(charge, mass, false);
            var decoyRankTable = _trainingParameters.GetRankTable(charge, mass, true);
            var prob = targetRankTable.GetRankProbability(rankNum, ionType).Found;
            var decoyProb = decoyRankTable.GetRankProbability(rankNum, ionType).Found;

            //if (prob == 0) Console.WriteLine("*** Target: {0} {1}", ionType.Name, rankNum);
            //if (decoyProb == 0) Console.WriteLine("*** Decoy: {0} {1}", ionType.Name, rankNum);
            if (Math.Abs(prob) < 1e-10 || Math.Abs(decoyProb) < 1e-10)
            {
                return 0;
            }

            return (Math.Log(prob) - Math.Log(decoyProb));
        }

        public IonType[] GetIonTypes(int charge, double mass)
        {
            return _trainingParameters.GetIonTypes(charge, mass);
        }

        private readonly TrainingParameters _trainingParameters;
    }
}
