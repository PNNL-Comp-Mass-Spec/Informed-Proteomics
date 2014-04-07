using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public abstract class IonFrequencyTable
    {
        private readonly Dictionary<string, IonProbability> _offsetFrequencies;

        protected IonFrequencyTable()
        {
            _offsetFrequencies = new Dictionary<string, IonProbability>();
        }

        public List<IonProbability> IonProbabilityTable
        {
            get { return _offsetFrequencies.Values.ToList(); }
        }

        public abstract void AddMatches(List<SpectrumMatch> matches);

        protected void AddIonProbability(IonProbability probability)
        {
            var name = probability.IonName;
            if (_offsetFrequencies.ContainsKey(name))
                _offsetFrequencies[name] += probability;
            else
                _offsetFrequencies.Add(name, probability);
        }
    }
}
