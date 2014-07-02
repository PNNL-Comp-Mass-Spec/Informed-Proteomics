using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDownViewer.Models
{
    public class PrSm: IComparable<PrSm>
    {
        public IcParameters Config { get; private set; }
        public int Scan { get; private set; }
        public ProductSpectrum Ms2Spectrum { get; private set; }
        public string Protein { get; private set; }
        public string Annotation { get; private set; }
        public string SequenceText { get; private set; }
        public List<string> SequenceLabel { get; private set; }
        public Sequence Sequence { get; private set; }
        public string Pre { get; private set; }
        public string Post { get; private set; }
        public List<Tuple<int, string>> Modifications { get; private set; }
        public string Composition { get; private set; }
        public string ProteinName { get; private set; }
        public string ProteinDesc { get; private set; }
        public int ProteinLength { get; private set; }
        public int Start { get; private set; }
        public int End { get; private set; }
        public int Charge { get; private set; }
        public double MostAbundantIsotopeMz { get; private set; }
        public double Mass { get; private set; }
        public int MatchedFragments { get; private set; }
        public double IsotopeCorrPrevMs1 { get; private set; }
        public double IsotopeCorrNextMs1 { get; private set; }
        public double CorrMostAbundantPlusOneIsoptope { get; private set; }
        public double ChargeCorrMinusOne { get; private set; }
        public double ChargeCorrPlusOne { get; private set; }
        public double QValue { get; private set; }
        public double PepQValue { get; private set; }

        public PrSm(string line, IcParameters config)
        {
            Config = config;
            _precursorXic = null;

            var parts = line.Split('\t');
            Scan = Convert.ToInt32(parts[0]);
            Ms2Spectrum = Config.Lcms.GetSpectrum(Scan) as ProductSpectrum;
            Pre = parts[1];
            Protein = parts[2];
            Post = parts[3];
            Annotation = (Pre + "." + Protein + "." + Post).Replace('-','_');
            SequenceLabel = new List<string>();
            SetModifications(parts[4]);
            Composition = parts[5];
            ProteinName = parts[6];
            ProteinDesc = parts[7].Split(';').FirstOrDefault();
            ProteinLength = Convert.ToInt32(parts[8]);
            Start = Convert.ToInt32(parts[9]);
            End = Convert.ToInt32(parts[10]);
            Charge = Convert.ToInt32(parts[11]);
            MostAbundantIsotopeMz = Convert.ToDouble(parts[12]);
            Mass = Convert.ToDouble(parts[13]);
            MatchedFragments = Convert.ToInt32(parts[14]);
            IsotopeCorrPrevMs1 = Convert.ToDouble(parts[15]);
            IsotopeCorrNextMs1 = Convert.ToDouble(parts[16]);
            CorrMostAbundantPlusOneIsoptope = Convert.ToDouble(parts[17]);
            ChargeCorrMinusOne = Convert.ToDouble(parts[18]);
            ChargeCorrPlusOne = Convert.ToDouble(parts[19]);
            QValue = Convert.ToDouble(parts[20]);
            PepQValue = Convert.ToDouble(parts[21]);
        }

        public PrSm() {}

        public List<Tuple<string, Peak[]>> GetIons(BaseIonType baseIon, int charge, List<Composition> fragments = null, List<int> indices = null)
        {
            var ionType = GetIonType(baseIon, NeutralLoss.NoLoss, charge);
            if (fragments == null) fragments = GetCompositions(ionType.BaseIonType.IsPrefix);
            if (indices != null && baseIon.IsPrefix) indices.Reverse();
            var peaks = new List<Tuple<string, Peak[]>>();
            for (int i = 0; i < fragments.Count; i++) 
            {
                var ion = ionType.GetIon(fragments[i]);
                var correlationScore = Ms2Spectrum.GetCorrScore(ion, Config.ProductIonTolerancePpm);
                if (correlationScore < CorrelationThreshold) continue;
                var isotopePeaks = Ms2Spectrum.GetAllIsotopePeaks(ion, Config.ProductIonTolerancePpm, RelativeIntensityThreshold);
                var baseIonTypeName = ionType.BaseIonType.Symbol;
                var index = i + 1;
                if (indices != null) index = indices[i];
                if (isotopePeaks != null && isotopePeaks.Length > 0)
                    peaks.Add(new Tuple<string, Peak[]>(baseIonTypeName + index + "-" + ionType.Charge + "+", isotopePeaks));
            }
            return peaks;
        }

        public List<Tuple<string, Xic>> GetFragmentIonXics(BaseIonType baseIon, int charge, List<Composition> fragments=null, List<int> indices=null)
        {
            var ionType = GetIonType(baseIon, NeutralLoss.NoLoss, charge);
            if (fragments == null) fragments = GetCompositions(ionType.BaseIonType.IsPrefix);
            var xics = new List<Tuple<string, Xic>>();
            if (indices != null && baseIon.IsPrefix) indices.Reverse();
            if (PrecursorXic == null || PrecursorXic.Count == 0) return xics;
            var precursorIon = Sequence.GetPrecursorIon(Charge);
            var precursorMz = precursorIon.GetMostAbundantIsotopeMz();

            for (int i = 0; i < fragments.Count; i++)
            {
                var ion = ionType.GetIon(fragments[i]);
                var correlationScore = Ms2Spectrum.GetCorrScore(ion, Config.ProductIonTolerancePpm);
                if (correlationScore < CorrelationThreshold) continue;
                var maxPeakMz = ion.GetMostAbundantIsotopeMz();
                var min = PrecursorXic.Min().ScanNum;
                var max = PrecursorXic.Max().ScanNum;
                var fragmentXic = Config.Lcms.GetExtractedFragmentIonChromatogram(maxPeakMz, Config.PrecursorTolerancePpm,
                                                                            precursorMz, min, max);
                var baseIonTypeName = ionType.BaseIonType.Symbol;
                var index = i + 1;
                if (indices != null) index = indices[i];
                xics.Add(new Tuple<string, Xic>(baseIonTypeName + index +"-"+ionType.Charge+"+", fragmentXic));
            }
            return xics;
        }

        public Xic PrecursorXic
        {
            get
            {
                if (_precursorXic == null)
                {
                    var precursorIon = Sequence.GetPrecursorIon(Charge);
                    var precursorMz = precursorIon.GetMostAbundantIsotopeMz();
                    _precursorXic = Config.Lcms.GetFullExtractedIonChromatogram(precursorMz, Config.PrecursorTolerancePpm);
                }
                return _precursorXic;
            }
        }

        public List<Tuple<string, Composition>> AnnotatedCompositions
        {
            get
            {
                var annotatedCompositions = new List<Tuple<string, Composition>>();
                for (int i = 0; i <= Sequence.Count; i++)
                {
                    var sequenceLabel = "";
                    if (i > 0) sequenceLabel = SequenceLabel[i - 1]; 
                    annotatedCompositions.Add(new Tuple<string, Composition>(sequenceLabel,
                                                                             Sequence.GetComposition(0, i)));
                }
                return annotatedCompositions;
            }
        }

        public int CompareTo(PrSm other)
        {
            return Scan.CompareTo(other.Scan);
        }

        private List<Composition> GetCompositions(bool prefix)
        {
            var compositions = new List<Composition>();
            for (int i = 1; i < Sequence.Count; i++)
            {
                compositions.Add(prefix
                    ? Sequence.GetComposition(0, i)
                    : Sequence.GetComposition(i, Sequence.Count));
            }
            if (!prefix) compositions.Reverse();
            return compositions;
        }

        private void SetModifications(string modifications)
        {
            Modifications = new List<Tuple<int, string>>();
            // Build Sequence AminoAcid list
            SequenceText = Protein;
            Sequence = new Sequence(SequenceText, new AminoAcidSet());
            foreach (var aa in Sequence) SequenceLabel.Add(aa.Residue.ToString(CultureInfo.InvariantCulture));

            // Parse modifications
            var mods = modifications.Split(',');
            if (mods.Length  <= 1) return;
            foreach (var modParts in mods.Select(mod => mod.Split(' ')))
            {
                if (modParts.Length < 0)    throw new FormatException("Invalid modification.");
                var modName = modParts[0];
                var modPos = Convert.ToInt32(modParts[1]);
                Modifications.Add(new Tuple<int, string>(modPos, modName));
            }
            
            // Add modifications to sequence
            Modifications.Sort(new CompareModByHighestPosition());   // sort in reverse order for insertion
            foreach (var mod in Modifications)
            {
                var pos = mod.Item1;
                if (pos > 0) pos--;
                var modLabel = "[" + mod.Item2 + "]";
                SequenceText = SequenceText.Insert(mod.Item1, modLabel);
                SequenceLabel[pos] += modLabel;
                var aa = Sequence[pos];
                var modaa = new ModifiedAminoAcid(aa, Modification.Get(mod.Item2));
                Sequence[pos] = modaa;
            }
        }

        private IonType GetIonType(BaseIonType baseIonType, NeutralLoss neutralLoss, int charge)
        {
            var chargeStr = charge.ToString(CultureInfo.InvariantCulture);
            if (charge == 1) chargeStr = "";
            var name = baseIonType.Symbol + chargeStr + neutralLoss.Name;
            return Config.IonTypeFactory.GetIonType(name);
        }

        private const double RelativeIntensityThreshold = 0.1;
        private const double CorrelationThreshold = 0.7;
        private Xic _precursorXic;
    }

    internal class CompareModByHighestPosition : IComparer<Tuple<int, string>>
    {
        public int Compare(Tuple<int, string> x, Tuple<int, string> y)
        {
            return (y.Item1.CompareTo(x.Item1));
        }
    }
}
