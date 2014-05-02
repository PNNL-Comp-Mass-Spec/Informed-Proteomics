using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
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

        public RankProbability(IEnumerable<IonType> ionTypes)
        {
            IonFrequencies = new Dictionary<IonType, double>();
            RankCount = 0;

            foreach (var ionType in ionTypes)
            {
                IonFrequencies.Add(ionType, 0);
            }
        }

        public void AddIons(Dictionary<IonType, double> ions)
        {
            foreach (var key in ions.Keys)
            {
                if (!IonFrequencies.ContainsKey(key))
                    IonFrequencies.Add(key, 0.0);
                IonFrequencies[key] += ions[key];
            }
        }
    }
}
