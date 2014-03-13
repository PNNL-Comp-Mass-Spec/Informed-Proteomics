using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

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

        public IEnumerable<IonProbability> GetCombinedChargeTable(IonTypeFactory ionTypeFactory)
        {
            var tempOff = new Dictionary<string, IonProbability>();

            foreach (var key in _offsetFrequencies.Keys)
            {
                if (_offsetFrequencies[key].Ion.Charge > 1)
                {
                    var reducedChargeName = key.Remove(1, 1);
                    if (!_offsetFrequencies.ContainsKey(reducedChargeName))
                        _offsetFrequencies.Add(reducedChargeName, new IonProbability(ionTypeFactory.GetIonType(reducedChargeName)));
                    tempOff[reducedChargeName] += _offsetFrequencies[key];
                }
                else
                {
                    tempOff.Add(key, _offsetFrequencies[key]);
                }
            }
            return tempOff.Values;
        }

        public void AddCleavageProbabilities(IEnumerable<SpectrumMatch> matches, IEnumerable<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold)
        {
            foreach (var match in matches)
            {
                foreach (var ionType in ionTypes)
                {
                    var prob = match.ContainsCleavageIons(ionType, tolerance, relativeIntensityThreshold);
                    AddIonProbability(prob);
                }
            }
        }

        public void AddPrecursorProbabilities(IEnumerable<SpectrumMatch> matches, IEnumerable<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold)
        {
            foreach (var match in matches)
            {
                foreach (var ionType in ionTypes)
                {
                    var conIon = match.ContainsPrecursorIon(ionType, tolerance, relativeIntensityThreshold);
                    int found = 0;
                    int total = 1;
                    if (conIon) found = 1;
                    var prob = new IonProbability(ionType, found, total); 
                    AddIonProbability(prob);
                }
            }
        }

        public void AddIonProbability(IonProbability probability)
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
