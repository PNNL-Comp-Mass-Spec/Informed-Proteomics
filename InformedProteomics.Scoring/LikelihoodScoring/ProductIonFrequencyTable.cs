using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class ProductIonFrequencyTable: IonFrequencyTable
    {
        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly double _relativeIntensityThreshold;
        private readonly bool _combineCharges;
        private readonly IonTypeFactory _ionTypeFactory;

        public ProductIonFrequencyTable(IEnumerable<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold, bool combineCharges=false)
        {
            _combineCharges = combineCharges;
            _ionTypes = ionTypes.ToArray();
            _tolerance = tolerance;
            _relativeIntensityThreshold = relativeIntensityThreshold;

            var charges = _ionTypes.Select(ionType => ionType.Charge).ToList();
            _ionTypeFactory = new IonTypeFactory(charges.Max());
        }

        private IonType ReducedChargeIonType(IonType ionType)
        {
            string name = ionType.Name;
            if (ionType.Charge > 1 && ionType.Charge < 10)
                name = name.Remove(1, 1);
            if (ionType.Charge >= 10)
                name = name.Remove(1, 2);
            return _ionTypeFactory.GetIonType(name);
        }

        public Probability<IonType>[] SelectIons(double probability)
        {
            var ionProbabilities = GetProbabilities();

            var selectedIons = ionProbabilities.Where(ionProbability => ionProbability.Prob >= probability).ToArray();
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
                    var ionTypeFound = new Dictionary<IonType, bool>();
                    foreach (var ionType in _ionTypes)
                    {
                        var it = ionType;
                        if (_combineCharges)
                            it = ReducedChargeIonType(ionType);
                        if (!ionTypeFound.ContainsKey(it))
                            ionTypeFound.Add(it, false);

                        var cleavagePoints = ionType.BaseIonType.IsPrefix ? prefixes : suffixes;

                        var ion = ionType.GetIon(cleavagePoints[i]);

                        if (spectrum.ContainsIon(ion, _tolerance, _relativeIntensityThreshold))
                        {
                            ionTypeFound[it] = true;
                        }
                    }
                    foreach (var key in ionTypeFound.Keys)
                    {
                        int found = 0;
                        const int total = 1;
                        if (ionTypeFound[key]) found = 1;
                        AddIonProbability(new Probability<IonType>(key, found, total));
                    }
                }
            }
        }

        public override IonType[] GetBinEdges()
        {
            return _ionTypes;
        }
    }
}
