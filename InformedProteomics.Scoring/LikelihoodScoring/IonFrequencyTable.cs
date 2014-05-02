using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public abstract class IonFrequencyTable: I1DProbabilityTable<IonType>
    {
        private readonly Dictionary<IonType, Probability<IonType>> _offsetFrequencies;

        protected IonFrequencyTable()
        {
            _offsetFrequencies = new Dictionary<IonType, Probability<IonType>>();
        }

        public abstract void AddMatches(List<SpectrumMatch> matches);

        protected void AddIonProbability(Probability<IonType> probability)
        {
            var name = probability.DataLabel;
            if (_offsetFrequencies.ContainsKey(name))
                _offsetFrequencies[name] += probability;
            else
                _offsetFrequencies.Add(name, probability);
        }

        public Probability<IonType>[] GetProbabilities()
        {
            return _offsetFrequencies.Values.ToArray();
        }

        public abstract IonType[] GetBinEdges();
    }
}
