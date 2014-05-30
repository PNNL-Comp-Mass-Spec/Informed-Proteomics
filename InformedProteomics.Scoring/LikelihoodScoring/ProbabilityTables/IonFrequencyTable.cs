using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public abstract class IonFrequencyTable
    {
        protected IonFrequencyTable()
        {
            _offsetFrequencies = new Dictionary<IonType, Probability<IonType>>();
        }

        public Probability<IonType>[] GetProbabilities()
        {
            return _offsetFrequencies.Values.ToArray();
        }

        public abstract void AddMatch(SpectrumMatch match);
        public abstract void AddMatches(List<SpectrumMatch> matches);
        public abstract void WriteToFile(StreamWriter file);
        public abstract IonType[] GetIonTypes();

        protected void AddIonProbability(Probability<IonType> probability)
        {
            var name = probability.Label;
            if (_offsetFrequencies.ContainsKey(name))
                _offsetFrequencies[name] += probability;
            else
                _offsetFrequencies.Add(name, probability);
        }
        private readonly Dictionary<IonType, Probability<IonType>> _offsetFrequencies;
    }
}
