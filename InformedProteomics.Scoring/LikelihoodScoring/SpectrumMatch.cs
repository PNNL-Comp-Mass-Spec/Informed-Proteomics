using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    class SpectrumMatch
    {
        private List<Composition> _prefixes;
        private List<Composition> _suffixes;
        private readonly Sequence _sequence;

        public string Peptide { get; private set; }
        public Spectrum Spectrum { get; private set; }

        public SpectrumMatch(string peptide, Spectrum spectrum, Sequence sequence)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            _sequence = sequence;
        }

        public SpectrumMatch(string protein, Spectrum spectrum)
        {
            Peptide = protein;
            Spectrum = spectrum;
            _sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(protein);
        }

        public Composition Comp
        {
            get { return _sequence.GetComposition(); }
        }

        public List<Composition> Prefixes
        {
            get
            {
                if (_prefixes == null)
                {
                    _prefixes = new List<Composition>();
                    for (int i = 1; i <= Peptide.Length; i++)
                    {
                        _prefixes.Add(_sequence.GetComposition(0, i));
                    }
                }
                return _prefixes;
            }
        }

        public List<Composition> Suffixes
        {
            get
            {
                if (_suffixes == null)
                {
                    _suffixes = new List<Composition>();
                    for (int i = 1; i <= Peptide.Length; i++)
                    {
                        _suffixes.Add(_sequence.GetComposition(Peptide.Length - i, Peptide.Length));
                    }
                }
                return _suffixes;
            }
        }

        public List<Ion> GetPrefixIons(IonType ionType)
        {
            return Prefixes.Select(ionType.GetIon).ToList();
        }

        public List<Ion> GetSuffixIons(IonType ionType)
        {
            return Suffixes.Select(ionType.GetIon).ToList();
        }

        public IonProbability ContainsCleavageIons(IonType ionType, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var ionTypeName = ionType.Name[0];
            var compositions = new List<Composition>();
            if (ionTypeName == 'a' || ionTypeName == 'b' || ionTypeName == 'c')
                compositions = Prefixes;
            else if (ionTypeName == 'x' || ionTypeName == 'y' || ionTypeName == 'z')
                compositions = Suffixes;

            var probability = new IonProbability(ionType);
            foreach (var composition in compositions)
            {
                var ion = ionType.GetIon(composition);
                probability.Total++;
                if (Spectrum.ContainsIon(ion, tolerance, relativeIntensityThreshold))
                    probability.Found++;
            }
            return probability;
        }
    }
}
