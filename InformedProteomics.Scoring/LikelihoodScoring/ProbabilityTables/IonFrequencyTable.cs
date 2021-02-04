using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class IonFrequencyTable
    {
        public IonFrequencyTable(IEnumerable<IonType> ionTypes, Tolerance tolerance, double relativeIntensityThreshold, bool combineCharges = false)
        {
            _ionFrequencies = new Dictionary<IonType, Probability<IonType>>();
            _combineCharges = combineCharges;
            _ionTypes = ionTypes.ToArray();
            _tolerance = tolerance;
            _relativeIntensityThreshold = relativeIntensityThreshold;

            var charges = _ionTypes.Select(ionType => ionType.Charge).ToList();
            _ionTypeFactory = new IonTypeFactory(charges.Max());
        }

        public Probability<IonType>[] GetProbabilities()
        {
            return _ionFrequencies.Values.ToArray();
        }

        /// <summary>
        /// Select ion types with probabilities at least as high as probabilityThreshold.
        /// </summary>
        /// <param name="probabilityThreshold"></param>
        /// <returns>Array of selected ion types.</returns>
        public IonType[] SelectIonTypes(double probabilityThreshold)
        {
            var ionProbabilities = GetProbabilities();
            return (from prob in ionProbabilities where prob.Prob >= probabilityThreshold select prob.Label).ToArray();
        }

        /// <summary>
        /// Add Peptide-Spectrum match to ion Probability table.
        /// </summary>
        /// <param name="match"></param>
        public void AddMatch(SpectrumMatch match)
        {
            var spectrum = match.Spectrum;

            var prefixes = match.Prefixes;
            var suffixes = match.Suffixes;

            for (var i = 0; i < prefixes.Count; i++)
            {
                var ionTypeFound = new Dictionary<IonType, bool>();
                foreach (var ionType in _ionTypes)
                {
                    var it = ionType;
                    if (_combineCharges)
                    {
                        it = ReducedChargeIonType(ionType);
                    }

                    if (!ionTypeFound.ContainsKey(it))
                    {
                        ionTypeFound.Add(it, false);
                    }

                    var cleavagePoints = ionType.BaseIonType.IsPrefix ? prefixes : suffixes;

                    var ion = ionType.GetIon(cleavagePoints[i]);

                    if (spectrum.ContainsIon(ion, _tolerance, _relativeIntensityThreshold))
                    {
                        ionTypeFound[it] = true;
                    }
                }
                foreach (var key in ionTypeFound.Keys)
                {
                    var found = 0;
                    const int total = 1;
                    if (ionTypeFound[key])
                    {
                        found = 1;
                    }

                    AddIonProbability(new Probability<IonType>(key, found, total));
                }
            }
        }

        /// <summary>
        /// Add list of Peptide-Spectrum matches to ion probability table.
        /// </summary>
        /// <param name="matches"></param>
        public void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                AddMatch(match);
            }
        }

        /// <summary>
        /// Write ion probability table to file.
        /// </summary>
        /// <param name="file">Open StreamWriter for file to write to.</param>
        public void WriteToFile(StreamWriter file)
        {
            var probabilities = GetProbabilities();
            foreach (var probability in probabilities)
            {
                file.WriteLine(probability.Label.Name + "\t" + probability.Prob);
            }
        }

        /// <summary>
        /// Get all possible ion types in ion probability table.
        /// </summary>
        /// <returns>Array of ion types.</returns>
        public IonType[] GetIonTypes()
        {
            return _ionTypes;
        }

        private void AddIonProbability(Probability<IonType> probability)
        {
            var name = probability.Label;
            if (_ionFrequencies.ContainsKey(name))
            {
                _ionFrequencies[name] += probability;
            }
            else
            {
                _ionFrequencies.Add(name, probability);
            }
        }

        private IonType ReducedChargeIonType(IonType ionType)
        {
            var name = ionType.Name;
            if (ionType.Charge > 1 && ionType.Charge < 10)
            {
                name = name.Remove(1, 1);
            }

            if (ionType.Charge >= 10)
            {
                name = name.Remove(1, 2);
            }

            return _ionTypeFactory.GetIonType(name);
        }

        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly double _relativeIntensityThreshold;
        private readonly bool _combineCharges;
        private readonly IonTypeFactory _ionTypeFactory;
        private readonly Dictionary<IonType, Probability<IonType>> _ionFrequencies;
    }
}
