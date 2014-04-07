using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class CleavageIonFrequencyTable: IonFrequencyTable
    {
        private readonly List<IonType> _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly double _relativeIntensityThreshold;
        private readonly bool _combineCharges;

        public CleavageIonFrequencyTable(List<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold, bool combineCharges=false)
        {
            _combineCharges = combineCharges;
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            _relativeIntensityThreshold = relativeIntensityThreshold;
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

        public IEnumerable<IonProbability> SelectIons(double probability)
        {
            var ionProbabilities = IonProbabilityTable;

            var selectedIons = ionProbabilities.Where(ionProbability => ionProbability.Prob >= probability).ToList();
            return selectedIons;
        }

        public override void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                var spectrum = match.Spectrum;

                var prefixes = match.Prefixes;
                var suffixes = match.Suffixes;

                for (int i = 0; i < prefixes.Count; i++)
                {
                    var ionTypeFound = new Dictionary<string, bool>();
                    foreach (var ionType in _ionTypes)
                    {
                        var name = ionType.Name;
                        if (_combineCharges)
                            name = ReducedChargeName(ionType);
                        if (!ionTypeFound.ContainsKey(name))
                            ionTypeFound.Add(name, false);

                        var cleavagePoints = ionType.BaseIonType.IsPrefix ? prefixes : suffixes;

                        var ion = ionType.GetIon(cleavagePoints[i]);

                        if (spectrum.ContainsIon(ion, _tolerance, _relativeIntensityThreshold))
                        {
                            ionTypeFound[name] = true;
                        }
                    }
                    foreach (var key in ionTypeFound.Keys)
                    {
                        int found = 0;
                        const int total = 1;
                        if (ionTypeFound[key]) found = 1;
                        AddIonProbability(new IonProbability(key, found, total));
                    }
                }
            }
        }
    }
}
