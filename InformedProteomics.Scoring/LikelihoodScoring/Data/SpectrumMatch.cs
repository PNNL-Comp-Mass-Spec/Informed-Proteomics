using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class SpectrumMatch
    {
        public Sequence Sequence { get; private set; }
        public string Peptide { get; private set; }
        public int PrecursorCharge { get; }
        public int ScanNum { get; }
        public Composition PrecursorComposition => Sequence.Composition + Composition.H2O;
        public bool Decoy { get; }

        public class MismatchException : Exception {
            public MismatchException() : base()
            {
            }

            public MismatchException(string message) : base(message)
            {
            }

            public MismatchException(string message, Exception innerException) : base(message, innerException)
            {
            }
        }

        public SpectrumMatch(Sequence sequence, Spectrum spectrum, int scanNum=0, int precursorCharge=1, bool decoy=false)
        {
            Peptide = string.Empty;
            foreach (var aa in sequence)
            {
                Peptide += aa.Residue;
            }

            _spectrum = spectrum;
            ScanNum = scanNum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            Sequence = sequence;
            if (decoy)
            {
                CreateDecoy();
            }
        }

        public SpectrumMatch(Sequence sequence, LazyLcMsRun lcms, int scanNum = 0, int precursorCharge = 1, bool decoy = false)
        {
            Peptide = string.Empty;
            foreach (var aa in sequence)
            {
                Peptide += aa.Residue;
            }

            _spectrum = null;
            _lcms = lcms;
            ScanNum = scanNum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            Sequence = sequence;
            if (decoy)
            {
                CreateDecoy();
            }
        }

        public SpectrumMatch(string peptide, DataFileFormat sequenceFormat,
                             LazyLcMsRun lcms, int scanNum = 0, int precursorCharge = 1,
                             bool decoy = false, string formula="")
        {
            Peptide = peptide;
            _spectrum = null;
            _lcms = lcms;
            ScanNum = scanNum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            var sequenceReader = new SequenceReader(sequenceFormat);
            Sequence = sequenceReader.GetSequence(peptide);
            if (decoy)
            {
                CreateDecoy();
            }
            else if (formula.Length > 0)
            {
                var composition = Composition.Parse(formula);
                if (!composition.Equals(Sequence.Composition + Composition.H2O))
                {
                    throw new MismatchException();
                }
            }
        }

        public SpectrumMatch(string peptide, DataFileFormat sequenceFormat,
                             Spectrum spectrum, int scanNum=0, int precursorCharge=1, bool decoy=false, string formula="")
        {
            Peptide = peptide;
            _spectrum = spectrum;
            _lcms = null;
            ScanNum = scanNum;
            PrecursorCharge = precursorCharge;
            Decoy = decoy;
            var sequenceReader = new SequenceReader(sequenceFormat);
            Sequence = sequenceReader.GetSequence(peptide);
            if (decoy)
            {
                CreateDecoy();
            }
            else if (formula.Length > 0)
            {
                var composition = Composition.Parse(formula);
                if (!composition.Equals(Sequence.Composition + Composition.H2O))
                {
                    throw new MismatchException();
                }
            }
        }

        public SpectrumMatch(SpectrumMatch match, bool decoy)
        {
            Peptide = match.Peptide;
            _spectrum = match._spectrum;
            _lcms = match._lcms;
            ScanNum = match.ScanNum;
            PrecursorCharge = match.PrecursorCharge;
            Decoy = decoy;
            Sequence = new Sequence(match.Sequence);
            if (decoy)
            {
                CreateDecoy();
            }
        }

        public Spectrum Spectrum => _spectrum ?? (_spectrum = _lcms.GetSpectrum(ScanNum));

        public Composition PeptideComposition => Sequence.Composition;

        public List<Composition> Prefixes
        {
            get
            {
                if (_prefixes == null)
                {
                    _prefixes = new List<Composition>();
                    for (var i = 1; i < Sequence.Count; i++)
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
                    for (var i = 1; i < Sequence.Count; i++)
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
            return compositions.ConvertAll(ionType.GetIon);
        }

        private void CreateDecoy()
        {
            Sequence.Reverse();
            var sequence = Sequence.Aggregate("", (current, aa) => current + aa.Residue);
            sequence = SimpleStringProcessing.Mutate(sequence, sequence.Length / 2);
            Peptide = sequence;
            Sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequence);
        }

        private Spectrum _spectrum;
        private readonly LazyLcMsRun _lcms;
        private List<Composition> _prefixes;
        private List<Composition> _suffixes;
    }
}
