using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class MassErrorScore
    {
        private readonly Dictionary<IonType, MassErrorTable> _trainingSet;
        private readonly Dictionary<IonType, MassErrorTable> _decoySet;

        public MassErrorScore(Dictionary<IonType, MassErrorTable> trainingSet, Dictionary<IonType, MassErrorTable> decoySet)
        {
            _trainingSet = trainingSet;
            _decoySet = decoySet;
        }

        public MassErrorScore(string fileName, int maxCharge = 2)
        {
            _trainingSet = new Dictionary<IonType, MassErrorTable>();
            _decoySet = new Dictionary<IonType, MassErrorTable>();
            ReadTrainingSetFromFile(fileName, maxCharge);
        }

        public double GetScore(IonType ionType, double value)
        {
            var targetProb = _trainingSet[ionType].GetMassErrorProbability(value).Found;
            if (targetProb.Equals(0.0)) targetProb = 1;
            var decoyProb = _decoySet[ionType].GetMassErrorProbability(value).Found;
            if (decoyProb.Equals(0.0)) decoyProb = 1;
            var score = (Math.Log(targetProb) - Math.Log(decoyProb));
            return score;
        }

        public IonType[] IonTypes
        {
            get { return _trainingSet.Keys.ToArray(); }
        }

        public void ReadTrainingSetFromFile(string fileName, int maxCharge, double offset=0.005)
        {
            var ionTypeFactory = new IonTypeFactory(maxCharge);
            var parser = new TsvFileParser(fileName);
            var data = parser.GetAllData();
            var ionTypes = (from key in data.Keys
                            where (key != "Error" && key != "Unexplained" &&
                                   key != "Total" && !key.Contains("De"))
                            select ionTypeFactory.GetIonType(key)).ToArray();
            var binEdges = data["Error"];
            var binEdgesD = binEdges.Select(binEdge => Convert.ToDouble(binEdge) - offset).ToArray();
            foreach (var ionType in ionTypes)
            {
                var targetData = data[ionType.Name].Select(Convert.ToDouble).ToList();
                var decoyData = data[ionType.Name + "-De"].Select(Convert.ToDouble).ToList();

                var targetMassErrorTable = new MassErrorTable(ionTypes, null, binEdgesD, targetData);
                var decoyMassErrorTable = new MassErrorTable(ionTypes, null, binEdgesD, decoyData);

                if (!_trainingSet.ContainsKey(ionType))
                    _trainingSet.Add(ionType, targetMassErrorTable);
                if (!_decoySet.ContainsKey(ionType))
                    _decoySet.Add(ionType, decoyMassErrorTable);
            }
        }
    }
}
