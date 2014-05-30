using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Config;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class SpectrumMatch
    {
        public Sequence Sequence { get; private set; }
        public string Peptide { get; private set; }
        public Spectrum Spectrum { get; private set; }
        public int PrecursorCharge { get; private set; }
        public Composition PrecursorComposition { get { return Sequence.Composition + Composition.H2O; }}
        public bool Decoy { get; private set; }
        public class MismatchException : Exception {}

        public SpectrumMatch(Sequence sequence, Spectrum spectrum, int precursorCharge=1, bool decoy=false)
        {
            foreach (var aa in sequence) Peptide += aa.Residue;
            Spectrum = spectrum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            Sequence = sequence;
            if (decoy) CreateDecoy();
        }

        public SpectrumMatch(string peptide, DataFileFormat sequenceFormat,
                             Spectrum spectrum, int precursorCharge=1, bool decoy=false)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            var sequenceReader = new SequenceReader(sequenceFormat);
            Sequence = sequenceReader.GetSequence(peptide);
            if (decoy) CreateDecoy();
        }

        public SpectrumMatch(string peptide, DataFileFormat sequenceFormat,
                             Spectrum spectrum, string formula, int precursorCharge=1, bool decoy=false)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            var sequenceReader = new SequenceReader(sequenceFormat);
            Sequence = sequenceReader.GetSequence(peptide);
            if (decoy) CreateDecoy();
            else
            {
                var composition = Composition.Parse(formula);
                if (!composition.Equals(Sequence.Composition + Composition.H2O))
                {
                    throw new MismatchException();
                }
            }
        }

        public Composition PeptideComposition
        {
            get { return Sequence.GetComposition(); }
        }

        public List<Composition> Prefixes
        {
            get
            {
                if (_prefixes == null)
                {
                    _prefixes = new List<Composition>();
                    for (int i = 1; i < Sequence.Count; i++)
                    {
                        _prefixes.Add(Sequence.GetComposition(0, i));
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
                    for (int i = 1; i < Sequence.Count; i++)
                    {
                        _suffixes.Add(Sequence.GetComposition(i, Sequence.Count));
                    }
                }
                return _suffixes;
            }
        }

        public List<Ion> GetCleavageIons(IonType ionType)
        {
            var compositions = ionType.BaseIonType.IsPrefix ? Prefixes : Suffixes;
            return compositions.Select(ionType.GetIon).ToList();
        }

        private void CreateDecoy()
        {
            Sequence.Reverse();
            var sequence = Sequence.Aggregate("", (current, aa) => current + aa.Residue);
            sequence = SimpleStringProcessing.Mutate(sequence, sequence.Length / 2);
            Sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequence);
        }

        private List<Composition> _prefixes;
        private List<Composition> _suffixes; 
    }
}
