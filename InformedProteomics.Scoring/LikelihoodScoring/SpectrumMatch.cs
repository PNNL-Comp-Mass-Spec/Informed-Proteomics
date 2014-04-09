using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumMatch
    {
        private List<Composition> _prefixes;
        private List<Composition> _suffixes;
        private readonly int _precursorCharge;

        public Sequence Sequence { get; private set; }
        public string Peptide { get; private set; }
        public Spectrum Spectrum { get; private set; }
        public int ScanNum { get; private set; }
        public int PrecursorCharge { get { return _precursorCharge;  } }

        public class MismatchException : Exception {}

        public SpectrumMatch(string peptide, Spectrum spectrum, int scanNum, int precursorCharge, Sequence sequence)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            ScanNum = scanNum;
            _precursorCharge = precursorCharge;
            Sequence = sequence;
        }

        public SpectrumMatch(string peptide, Spectrum spectrum, int scanNum, int precursorCharge)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            ScanNum = scanNum;
            _precursorCharge = precursorCharge;
            Sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(peptide);
        }

        public SpectrumMatch(string peptide, Spectrum spectrum, int scanNum, int precursorCharge, string formula)
        {
            Peptide = peptide;
            Spectrum = spectrum;
            ScanNum = scanNum;
            _precursorCharge = precursorCharge;
            Sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(peptide);

            var composition = Composition.Parse(formula);
            if (!composition.Equals(Sequence.Composition + Composition.H2O))
            {
                throw new MismatchException();
            }
        }

        public Composition PeptideComposition
        {
            get { return Sequence.GetComposition(); }
        }

        public string GetCleanPeptide()
        {
            var aaset = new AminoAcidSet();
            var cleanPeptide = "";
            foreach (var aminoacid in Peptide)
            {
                if (aaset.GetAminoAcid(aminoacid) != null)
                    cleanPeptide += aminoacid;
            }
            return cleanPeptide;
        }

        public string GetPeptidePrefix()
        {
            var cleanPeptide = GetCleanPeptide();
            return cleanPeptide.Substring(0, cleanPeptide.Length - 1);
        }

        public string GetPeptideSuffix()
        {
            var cleanPeptide = GetCleanPeptide();
            return cleanPeptide.Substring(1, cleanPeptide.Length);
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
                        _suffixes.Add(Sequence.GetComposition(Peptide.Length - i, Peptide.Length));
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

        public List<Ion> GetCleavageIons(IonType ionType)
        {
            var compositions = ionType.BaseIonType.IsPrefix ? Prefixes : Suffixes;
            return compositions.Select(ionType.GetIon).ToList();
        }

        public IonProbability ContainsCleavageIons(IonType ionType, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var ions = GetCleavageIons(ionType);

            var probability = new IonProbability(ionType.Name);
            foreach (var ion in ions)
            {
                probability.Total++;
                if (Spectrum.ContainsIon(ion, tolerance, relativeIntensityThreshold))
                    probability.Found++;
            }
            return probability;
        }

        public bool ContainsPrecursorIon(IonType ionType, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var ion = ionType.GetIon(PeptideComposition);
            return Spectrum.ContainsIon(ion, tolerance, relativeIntensityThreshold);
        }
    }
}
