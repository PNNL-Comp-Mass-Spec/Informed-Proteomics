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

        public List<IonProbability> IonProbabilityTable
        {
            get { return _offsetFrequencies.Values.ToList(); }
        }

        private string ReducedChargeName(IonType ionType)
        {
            string name = ionType.Name;
            if (ionType.Charge > 1 && ionType.Charge < 10)
                name = name.Remove(1, 1);
            if (ionType.Charge >= 10)
                name = name.Remove(1, 2);
            return name;
        }

        public void AddCleavageProbabilities(List<SpectrumMatch> matches, List<IonType> ionTypes,
            Tolerance tolerance, double relativeIntensityThreshold, bool combineCharges=false)
        {
            foreach (var match in matches)
            {
                var spectrum = match.Spectrum;

                var prefixes = match.Prefixes;
                var suffixes = match.Suffixes;

                for (int i = 0; i < prefixes.Count; i++)
                {
                    var ionTypeFound = new Dictionary<string, bool>();
                    foreach (var ionType in ionTypes)
                    {
                        var name = ionType.Name;
                        if (combineCharges)
                            name = ReducedChargeName(ionType);
                        if (!ionTypeFound.ContainsKey(name))
                            ionTypeFound.Add(name, false);

                        var cleavagePoints = ionType.BaseIonType.IsPrefix ? prefixes : suffixes;

                        var ion = ionType.GetIon(cleavagePoints[i]);

                        if (spectrum.ContainsIon(ion, tolerance, relativeIntensityThreshold))
                        {
                            ionTypeFound[name] = true;
                        }
/*                        if (spectrum.FindPeak(ion.GetMonoIsotopicMz(), tolerance) != null)
                        {
                            ionTypeFound[name] = true;
                        }*/
                    }
                    foreach (var key in ionTypeFound.Keys)
                    {
                        int found = 0;
                        int total = 1;
                        if (ionTypeFound[key]) found = 1;
                        AddIonProbability(new IonProbability(key, found, total));
                    }
                }
            }
        }

        public void AddPrecursorProbabilities(List<SpectrumMatch> matches, List<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold)
        {
            foreach (var match in matches)
            {
                foreach (var ionType in ionTypes)
                {
                    var conIon = match.ContainsPrecursorIon(ionType, tolerance, relativeIntensityThreshold);
                    int found = 0;
                    int total = 1;
                    if (conIon) found = 1;
                    var prob = new IonProbability(ionType.Name, found, total); 
                    AddIonProbability(prob);
                }
            }
        }

        public void AddIonProbability(IonProbability probability)
        {
            var name = probability.IonName;
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
