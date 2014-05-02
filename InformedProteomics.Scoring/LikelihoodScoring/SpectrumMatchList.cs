using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumMatchList: List<SpectrumMatch>
    {
        private const string QValueHeader = "Qvalue";
        private const string EvalueHeader = "E-value";
        private const string TopDownPeptideHeader = "Annotation";
        private const double QValueThreshold = 0.01;
        private const double EValueThreshold = 0.1;

        private const string BottomUpPeptideHeader = "Peptide";
        private const string PepQValueHeader = "PepQValue";
        private const string FormulaHeader = "Formula";
        private const double PepQValueThreshold = 0.01;

        private const string ScanHeader = "ScanNum";
        private const string PrecursorChargeHeader = "Charge";

        public int MaxCharge { get; private set; }
        public bool Decoy { get; private set; }
        public ActivationMethod Act { get; private set; }

        public SpectrumMatchList(LcMsRun lcms, TsvFileParser tsvFile, ActivationMethod act, bool useDecoy=false, int maxCharge=0)
        {
            Decoy = useDecoy;
            Act = act;
            MaxCharge = maxCharge;
            AddMatchesFromFile(lcms, tsvFile);
        }

        public SpectrumMatchList(ActivationMethod act, bool useDecoy=false, int maxCharge=0)
        {
            Decoy = useDecoy;
            MaxCharge = maxCharge;
            Act = act;
        }

        public void AddMatchesFromFile(LcMsRun lcms, TsvFileParser tsvFile)
        {
            var precursorCharges = tsvFile.GetData(PrecursorChargeHeader);
            var scans = tsvFile.GetData(ScanHeader);

            var peptides = tsvFile.GetData(TopDownPeptideHeader);
            if (peptides != null)
            {
                var peptideSet = new HashSet<string>();
                var filterThreshold = QValueThreshold;
                var filterValues = tsvFile.GetData(QValueHeader);
                if (filterValues == null)
                {
                    filterValues = tsvFile.GetData(EvalueHeader);
                    filterThreshold = EValueThreshold;
                }

                var aset = new AminoAcidSet();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(filterValues[i]) > filterThreshold || peptideSet.Contains(peptides[i])) continue;
                    peptideSet.Add(peptides[i]);
                    var scanNum = Convert.ToInt32(scans[i]);
                    var spectrum = lcms.GetSpectrum(scanNum);
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != Act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    if (MaxCharge > 0 && precursorCharge > MaxCharge) continue;
                    Add(new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge, new Sequence(peptides[i], aset), Decoy));
                }
            }
            else
            {
                peptides = tsvFile.GetData(BottomUpPeptideHeader);
                if (scans == null)
                    throw new FormatException();

                var pepQValues = tsvFile.GetData(PepQValueHeader);
                var formulas = tsvFile.GetData(FormulaHeader);

                var peptideSet = new HashSet<string>();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(pepQValues[i]) > PepQValueThreshold || peptideSet.Contains(peptides[i])) continue;
                    peptideSet.Add(peptides[i]);
                    var scanNum = Convert.ToInt32(scans[i]);
                    var spectrum = lcms.GetSpectrum(scanNum);
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != Act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    if (MaxCharge > 0 && precursorCharge > MaxCharge) continue;
                    Add((formulas != null && formulas[i] != null)
                        ? new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge, formulas[i], Decoy)
                        : new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge, Decoy));
                }
            }
        }

        public void AddMatch(SpectrumMatch match)
        {
            var newMatch = new SpectrumMatch(match.Peptide, match.Spectrum, match.ScanNum,
                                 match.PrecursorCharge, Decoy);
            Add(newMatch);
        }

        public void AddMatches(IEnumerable<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                AddMatch(match);
            }
        }

        public void FilterSpectra(double windowWidth=100, int retentionCount=6)
        {
            var filteredList = new SpectrumMatchList(Act);
            foreach (var match in this)
            {
                var spectrum = SpectrumFilter.GetFilteredSpectrum(match.Spectrum, windowWidth, retentionCount);
                filteredList.Add(new SpectrumMatch(match.Peptide, spectrum, match.ScanNum, match.PrecursorCharge));
            }
            Clear();
            AddRange(filteredList);
        }

        public SpectrumMatch GetScan(int scanNum)
        {
            return (from i in this
                    where i.ScanNum == scanNum
                    select i).FirstOrDefault();
        }

        public SpectrumMatchList GetCharge(int charge)
        {
            var chargeMatchList = new SpectrumMatchList(Act);
            chargeMatchList.AddRange(from i in this
                    where i.PrecursorCharge == charge
                    select i);
            return chargeMatchList;
        }
    }
}
