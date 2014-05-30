using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class RankProbability
    {
        public Dictionary<IonType, double> IonFrequencies { get; private set; }
        public int RankCount { get; set; }

        public Dictionary<IonType, Probability<string>> IonProbabilities
        {
            get
            {
                var probDict = IonFrequencies.Keys.ToDictionary(ionType => ionType, ionType => new Probability<string>(ionType.Name, IonFrequencies[ionType], RankCount));
                return probDict;
            }
        }

        public RankProbability()
        {
            IonFrequencies = new Dictionary<IonType, double>();
            RankCount = 0;
        }

        public RankProbability(IEnumerable<IonType> ionTypes)
        {
            IonFrequencies = new Dictionary<IonType, double>();
            foreach (var ionType in ionTypes) IonFrequencies.Add(ionType, 0.0);
            RankCount = 0;
        }

        public void AddIon(IonType ionType, double found)
        {
            if (!IonFrequencies.ContainsKey(ionType))
                IonFrequencies.Add(ionType, 0.0);
            IonFrequencies[ionType] += found;
        }

        public void AddIons(Dictionary<IonType, double> ions)
        {
            foreach (var key in ions.Keys)
            {
                AddIon(key, ions[key]);
            }
        }
    }
}
