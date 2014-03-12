using System.Collections.Generic;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class OffsetFrequencyTable
    {
        private readonly Dictionary<string, IonProbability> _offsetFrequencies;

        public OffsetFrequencyTable()
        {
            _offsetFrequencies = new Dictionary<string, IonProbability>();
        }

        public IEnumerable<IonProbability> IonProbabilityTable
        {
            get { return _offsetFrequencies.Values; }
        }

        public IEnumerable<IonProbability> CombinedChargeTable
        {
            get
            {
                var tempOff = new Dictionary<string, IonProbability>();

                foreach (var key in _offsetFrequencies.Keys)
                {
                    if (_offsetFrequencies[key].Ion.Charge > 1)
                    {
                        var reducedChargeName = key.Remove(1, 1);
                        if (!_offsetFrequencies.ContainsKey(reducedChargeName))
                            _offsetFrequencies.Add(reducedChargeName, new IonProbability(_offsetFrequencies[key].Ion));
                        tempOff[reducedChargeName] += _offsetFrequencies[key];
                    }
                    else
                    {
                        tempOff.Add(key, _offsetFrequencies[key]);
                    }
                }
                return tempOff.Values;
            }
        }

        public void AddOffsetFrequency(IonProbability probability)
        {
            var name = probability.Ion.Name;
            if (_offsetFrequencies.ContainsKey(name))
            {
                _offsetFrequencies[name].Found += probability.Found;
                _offsetFrequencies[name].Total += probability.Total;
            }
            else
            {
                _offsetFrequencies.Add(name, probability);
            }
        }


    }
}
