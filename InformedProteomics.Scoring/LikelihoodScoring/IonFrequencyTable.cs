using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public abstract class IonFrequencyTable
    {
        private readonly Dictionary<string, Probability<string>> _offsetFrequencies;

        protected IonFrequencyTable()
        {
            _offsetFrequencies = new Dictionary<string, Probability<string>>();
        }

        public List<Probability<string>> IonProbabilityTable
        {
            get { return _offsetFrequencies.Values.ToList(); }
        }

        public abstract void AddMatches(List<SpectrumMatch> matches);

        protected void AddIonProbability(Probability<string> probability)
        {
            var name = probability.DataLabel;
            if (_offsetFrequencies.ContainsKey(name))
                _offsetFrequencies[name] += probability;
            else
                _offsetFrequencies.Add(name, probability);
        }
    }
}
